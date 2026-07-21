// Dependency-free print remeshing: generalized-winding-number signed distance
// sampled on a regular grid, contoured with the marching-tetrahedra extractor
// (examples/marching_cubes.h). Unlike the CGAL Alpha Wrap backend this tracks
// the input surface exactly (no offset dilation) wherever the soup encloses
// volume, at the cost of grid-limited resolution. Open sheets have no interior
// under a winding classification, so ambiguous regions (winding near 1/2) fall
// back to a thin offset shell around the surface — closed geometry is left
// untouched, sheets become printable slabs about one voxel thick.
//
// Everything here is self-contained: one triangle BVH provides both exact
// closest-point queries (the distance term, also reused for PBR projection)
// and a Barnes-Hut evaluation of the winding number [Barill et al. 2018,
// first-order far field, exact Van Oosterom-Strackee solid angles up close].

#include "print_remesh.h"

#include "marching_cubes.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>

namespace t2print {

namespace {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

inline Vec3 operator-(const Vec3 & a, const Vec3 & b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator+(const Vec3 & a, const Vec3 & b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator*(double s, const Vec3 & a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(const Vec3 & a, const Vec3 & b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3 cross(const Vec3 & a, const Vec3 & b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
inline double norm(const Vec3 & a) { return std::sqrt(dot(a, a)); }

// Ericson's exact closest point on a triangle (Real-Time Collision Detection,
// 5.1.5). Also reports the clamped barycentric coordinates for interpolation.
Vec3 closest_on_triangle(const Vec3 & p, const Vec3 & a, const Vec3 & b, const Vec3 & c,
                         double bary[3]) {
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) { bary[0]=1; bary[1]=0; bary[2]=0; return a; }
    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) { bary[0]=0; bary[1]=1; bary[2]=0; return b; }
    const double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        bary[0]=1-v; bary[1]=v; bary[2]=0; return a + v*ab;
    }
    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) { bary[0]=0; bary[1]=0; bary[2]=1; return c; }
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        bary[0]=1-w; bary[1]=0; bary[2]=w; return a + w*ac;
    }
    const double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bary[0]=0; bary[1]=1-w; bary[2]=w; return b + w*(c - b);
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    bary[0]=1-v-w; bary[1]=v; bary[2]=w;
    return a + v*ab + w*ac;
}

struct BVH {
    struct Node {
        double bmin[3], bmax[3];
        int32_t left = -1, right = -1;  // internal children, or
        int32_t start = 0, count = 0;   // leaf range into `order`
        Vec3 area_normal;               // sum of triangle area-weighted normals
        Vec3 centroid;                  // area-weighted centroid (winding dipole)
        double area = 0.0;
        double radius = 0.0;            // covers every triangle vertex from centroid
    };

    std::vector<Vec3>    a, b, c;       // triangle vertices, in `order` order
    std::vector<int32_t> order;         // node ranges -> original triangle index
    std::vector<Node>    nodes;

    static constexpr int    LEAF = 8;
    static constexpr double BETA2 = 4.0;   // Barnes-Hut far-field threshold, beta = 2

    void build(const std::vector<float> & verts, const std::vector<int32_t> & tris) {
        const size_t nt = tris.size() / 3;
        std::vector<Vec3> ta(nt), tb(nt), tc(nt), tm(nt);
        order.reserve(nt);
        for (size_t t = 0; t < nt; ++t) {
            const float * pa = verts.data() + (size_t) tris[3*t]   * 3;
            const float * pb = verts.data() + (size_t) tris[3*t+1] * 3;
            const float * pc = verts.data() + (size_t) tris[3*t+2] * 3;
            ta[t] = {pa[0], pa[1], pa[2]};
            tb[t] = {pb[0], pb[1], pb[2]};
            tc[t] = {pc[0], pc[1], pc[2]};
            tm[t] = (1.0/3.0) * (ta[t] + tb[t] + tc[t]);
            order.push_back((int32_t) t);
        }
        nodes.reserve(2 * nt / LEAF + 4);
        build_node(ta, tb, tc, tm, 0, (int32_t) nt);
        a.resize(nt); b.resize(nt); c.resize(nt);
        for (size_t i = 0; i < nt; ++i) {
            a[i] = ta[order[i]]; b[i] = tb[order[i]]; c[i] = tc[order[i]];
        }
    }

    int32_t build_node(const std::vector<Vec3> & ta, const std::vector<Vec3> & tb,
                       const std::vector<Vec3> & tc, const std::vector<Vec3> & tm,
                       int32_t start, int32_t count) {
        const int32_t ni = (int32_t) nodes.size();
        nodes.emplace_back();
        {
            Node & n = nodes[ni];
            n.bmin[0] = n.bmin[1] = n.bmin[2] = 1e300;
            n.bmax[0] = n.bmax[1] = n.bmax[2] = -1e300;
            for (int32_t i = start; i < start + count; ++i) {
                const int32_t t = order[i];
                for (const Vec3 * v : {&ta[t], &tb[t], &tc[t]}) {
                    const double p[3] = {v->x, v->y, v->z};
                    for (int k = 0; k < 3; ++k) {
                        n.bmin[k] = std::min(n.bmin[k], p[k]);
                        n.bmax[k] = std::max(n.bmax[k], p[k]);
                    }
                }
            }
        }
        if (count <= LEAF) {
            Node & n = nodes[ni];
            n.start = start; n.count = count;
            Vec3 wc{};
            for (int32_t i = start; i < start + count; ++i) {
                const int32_t t = order[i];
                const Vec3 an = 0.5 * cross(tb[t] - ta[t], tc[t] - ta[t]);
                const double area = norm(an);
                n.area_normal = n.area_normal + an;
                n.area += area;
                wc = wc + area * tm[t];
            }
            n.centroid = n.area > 0.0 ? (1.0 / n.area) * wc
                                      : Vec3{0.5*(n.bmin[0]+n.bmax[0]),
                                             0.5*(n.bmin[1]+n.bmax[1]),
                                             0.5*(n.bmin[2]+n.bmax[2])};
            for (int32_t i = start; i < start + count; ++i) {
                const int32_t t = order[i];
                for (const Vec3 * v : {&ta[t], &tb[t], &tc[t]})
                    n.radius = std::max(n.radius, norm(*v - n.centroid));
            }
            return ni;
        }
        int axis = 0;
        {
            const Node & n = nodes[ni];
            const double ext[3] = {n.bmax[0]-n.bmin[0], n.bmax[1]-n.bmin[1], n.bmax[2]-n.bmin[2]};
            if (ext[1] > ext[axis]) axis = 1;
            if (ext[2] > ext[axis]) axis = 2;
        }
        const int32_t mid = start + count / 2;
        auto key = [&](int32_t t) {
            return axis == 0 ? tm[t].x : axis == 1 ? tm[t].y : tm[t].z;
        };
        std::nth_element(order.begin() + start, order.begin() + mid,
                         order.begin() + start + count,
                         [&](int32_t l, int32_t r) { return key(l) < key(r); });
        const int32_t li = build_node(ta, tb, tc, tm, start, mid - start);
        const int32_t ri = build_node(ta, tb, tc, tm, mid, start + count - mid);
        Node & n = nodes[ni];
        n.left = li; n.right = ri;
        const Node & l = nodes[li], & r = nodes[ri];
        n.area_normal = l.area_normal + r.area_normal;
        n.area = l.area + r.area;
        n.centroid = n.area > 0.0
            ? (1.0 / n.area) * (l.area * l.centroid + r.area * r.centroid)
            : Vec3{0.5*(n.bmin[0]+n.bmax[0]), 0.5*(n.bmin[1]+n.bmax[1]), 0.5*(n.bmin[2]+n.bmax[2])};
        n.radius = std::max(norm(l.centroid - n.centroid) + l.radius,
                            norm(r.centroid - n.centroid) + r.radius);
        return ni;
    }

    double box_dist2(const Node & n, const Vec3 & p) const {
        double d2 = 0.0;
        const double q[3] = {p.x, p.y, p.z};
        for (int k = 0; k < 3; ++k) {
            const double d = std::max({n.bmin[k] - q[k], 0.0, q[k] - n.bmax[k]});
            d2 += d * d;
        }
        return d2;
    }

    // Closest point on the soup; returns squared distance and fills the winning
    // triangle (original index) and its clamped barycentric coordinates.
    double closest(const Vec3 & p, int32_t & tri, double bary[3]) const {
        double best2 = 1e300;
        tri = -1;
        int32_t stack[128];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node & n = nodes[stack[--sp]];
            if (box_dist2(n, p) >= best2) continue;
            if (n.left < 0) {
                for (int32_t i = n.start; i < n.start + n.count; ++i) {
                    double bc[3];
                    const Vec3 q = closest_on_triangle(p, a[i], b[i], c[i], bc);
                    const double d2 = dot(q - p, q - p);
                    if (d2 < best2) {
                        best2 = d2; tri = order[i];
                        bary[0] = bc[0]; bary[1] = bc[1]; bary[2] = bc[2];
                    }
                }
                continue;
            }
            // Descend the nearer child first so the farther one usually prunes.
            const double dl = box_dist2(nodes[n.left], p), dr = box_dist2(nodes[n.right], p);
            if (dl < dr) { stack[sp++] = n.right; stack[sp++] = n.left; }
            else         { stack[sp++] = n.left;  stack[sp++] = n.right; }
        }
        return best2;
    }

    // Generalized winding number: exact signed solid angles near the query,
    // single-dipole approximation for well-separated clusters.
    double winding(const Vec3 & p) const {
        constexpr double INV4PI = 1.0 / (4.0 * 3.14159265358979323846);
        double w = 0.0;
        int32_t stack[128];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node & n = nodes[stack[--sp]];
            const Vec3 d = n.centroid - p;
            const double d2 = dot(d, d);
            if (d2 > BETA2 * n.radius * n.radius && d2 > 0.0) {
                w += INV4PI * dot(n.area_normal, d) / (d2 * std::sqrt(d2));
                continue;
            }
            if (n.left >= 0) {
                stack[sp++] = n.left;
                stack[sp++] = n.right;
                continue;
            }
            for (int32_t i = n.start; i < n.start + n.count; ++i) {
                // Van Oosterom & Strackee's signed solid angle of a triangle.
                const Vec3 va = a[i] - p, vb = b[i] - p, vc = c[i] - p;
                const double la = norm(va), lb = norm(vb), lc = norm(vc);
                const double det = dot(va, cross(vb, vc));
                const double den = la*lb*lc + dot(va, vb)*lc + dot(vb, vc)*la + dot(vc, va)*lb;
                w += 2.0 * INV4PI * std::atan2(det, den);
            }
        }
        return w;
    }
};

void parallel_for(size_t total, size_t chunk, const std::function<void(size_t, size_t)> & fn) {
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    const unsigned workers = (unsigned) std::min<size_t>(std::min(16u, hw), (total + chunk - 1) / chunk);
    if (workers <= 1) { fn(0, total); return; }
    std::atomic<size_t> next{0};
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for (unsigned i = 0; i < workers; ++i)
        threads.emplace_back([&]() {
            for (;;) {
                const size_t begin = next.fetch_add(chunk);
                if (begin >= total) break;
                fn(begin, std::min(total, begin + chunk));
            }
        });
    for (auto & t : threads) t.join();
}

// Drop output islands below min_frac of the triangles (the winding field can
// leave speck "bubbles" where a messy soup makes the classification ambiguous),
// mirroring mesh_export's filter_components. Always keeps the largest island.
void drop_small_islands(std::vector<float> & verts, std::vector<float> & normals,
                        std::vector<int32_t> & tris, float min_frac) {
    const int32_t nv = (int32_t)(verts.size() / 3), nt = (int32_t)(tris.size() / 3);
    if (nv == 0 || nt == 0) return;
    std::vector<int32_t> parent((size_t) nv);
    for (int32_t i = 0; i < nv; ++i) parent[i] = i;
    auto find = [&](int32_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (int32_t t = 0; t < nt; ++t) {
        const int32_t r = find(tris[3*t]);
        parent[find(tris[3*t+1])] = r;
        parent[find(tris[3*t+2])] = r;
    }
    std::vector<int32_t> tris_in((size_t) nv, 0);
    int32_t largest = -1;
    for (int32_t t = 0; t < nt; ++t) {
        const int32_t r = find(tris[3*t]);
        if (++tris_in[r] > (largest < 0 ? 0 : tris_in[largest])) largest = r;
    }
    const int32_t keep_min = std::max(1, (int32_t)(nt * min_frac));
    std::vector<int32_t> vremap((size_t) nv, -1);
    std::vector<float> nvts, nnorm;
    std::vector<int32_t> nfs;
    for (int32_t t = 0; t < nt; ++t) {
        const int32_t r = find(tris[3*t]);
        if (r != largest && tris_in[r] < keep_min) continue;
        for (int k = 0; k < 3; ++k) {
            const int32_t vi = tris[3*t+k];
            if (vremap[vi] < 0) {
                vremap[vi] = (int32_t)(nvts.size() / 3);
                nvts.insert(nvts.end(), verts.begin() + (size_t) vi*3, verts.begin() + (size_t) vi*3 + 3);
                nnorm.insert(nnorm.end(), normals.begin() + (size_t) vi*3, normals.begin() + (size_t) vi*3 + 3);
            }
            nfs.push_back(vremap[vi]);
        }
    }
    verts.swap(nvts); normals.swap(nnorm); tris.swap(nfs);
}

} // namespace

bool winding_remesh(const std::vector<float> & source_verts,
                    const std::vector<int32_t> & source_tris,
                    int grid_res,
                    float shell_ratio,
                    std::vector<float> & out_verts,
                    std::vector<float> & out_normals,
                    std::vector<int32_t> & out_tris,
                    std::string & err) {
    out_verts.clear(); out_normals.clear(); out_tris.clear();
    if (source_verts.size() < 9 || source_tris.size() < 3) {
        err = "empty mesh";
        return false;
    }
    grid_res = std::max(32, std::min(1024, grid_res));

    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (size_t i = 0; i < source_verts.size() / 3; ++i) {
        const float * p = source_verts.data() + i * 3;
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
            err = "mesh contains a non-finite vertex";
            return false;
        }
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], (double) p[k]);
            hi[k] = std::max(hi[k], (double) p[k]);
        }
    }
    const size_t nv = source_verts.size() / 3;
    std::vector<int32_t> tris;
    tris.reserve(source_tris.size());
    for (size_t t = 0; t < source_tris.size() / 3; ++t) {
        const int32_t a = source_tris[3*t], b = source_tris[3*t+1], c = source_tris[3*t+2];
        if (a < 0 || b < 0 || c < 0 || (size_t) a >= nv || (size_t) b >= nv || (size_t) c >= nv) {
            err = "triangle index out of range";
            return false;
        }
        if (a == b || b == c || a == c) continue;
        tris.insert(tris.end(), {a, b, c});
    }
    if (tris.empty()) { err = "mesh has no valid triangles"; return false; }

    const double ext[3] = {hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]};
    const double diagonal = std::sqrt(ext[0]*ext[0] + ext[1]*ext[1] + ext[2]*ext[2]);
    if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
        err = "mesh has an empty bounding box";
        return false;
    }

    // Grid: `grid_res` cells along the longest axis, three-plus padding cells
    // per side so the surface (and any shell around it) never touches the
    // border. The extra irrational origin shift keeps input geometry (often
    // axis-aligned) off the corner lattice: a surface running exactly through
    // grid corners makes distance, winding and field gradient all degenerate
    // there, which the contouring's tie-breaking cannot always untangle.
    const double h = std::max({ext[0], ext[1], ext[2]}) / grid_res;
    double glo[3];
    int dims[3];
    for (int k = 0; k < 3; ++k) {
        glo[k] = lo[k] - (3.0 + 0.6180339887) * h;
        dims[k] = (int) std::ceil(ext[k] / h) + 9;   // cells + ~7.6 padding, +1 corner
    }
    const size_t total = (size_t) dims[0] * dims[1] * dims[2];
    if (total > 400u * 1000 * 1000) {
        err = "print remesh grid too large (lower T2_PRINT_RES)";
        return false;
    }
    // Open sheets survive as slabs of at least one voxel; anything thinner
    // cannot be represented by the grid (and would not be printable anyway).
    const bool  shell_on = shell_ratio >= 0.0f;
    const double shell_t = std::max((double) shell_ratio * diagonal, 0.9 * h);
    // Exact-evaluation band around the surface (see the coarse/fine fill
    // below); also the distance over which the shell term tapers out.
    const double band = std::max(7.0 * h, 2.0 * shell_t + 3.0 * h);

    BVH bvh;
    bvh.build(source_verts, tris);

    // Exact winding-signed distance with an open-sheet fallback: when a point
    // near the surface is classified outside AND so is its mirror through the
    // closest surface point, the soup encloses volume on neither side — an open
    // sheet (or its rim), which the plain sign would erase. There, and only
    // there, fall back to an offset shell. Near closed surfaces the mirror
    // lands inside, so their zero level is never dilated; at sheet rims the
    // mirror stays outside in every direction, so the shell wraps the rim
    // instead of pinching to a non-manifold knife edge.
    auto sample = [&](const Vec3 & p) -> float {
        int32_t tri; double bary[3];
        const double d = std::sqrt(bvh.closest(p, tri, bary));
        const double w = std::fabs(bvh.winding(p));
        double s = w > 0.5 ? -d : d;
        if (shell_on && w <= 0.5 && d < band) {
            const int32_t * ids = tris.data() + (size_t) tri * 3;
            Vec3 c{};
            for (int k = 0; k < 3; ++k) {
                const float * v = source_verts.data() + (size_t) ids[k] * 3;
                c = c + bary[k] * Vec3{v[0], v[1], v[2]};
            }
            const Vec3 mirror = 2.0 * c - p;
            if (std::fabs(bvh.winding(mirror)) <= 0.5) {
                // Taper the shell to nothing at the band edge so the field
                // stays continuous where the fill switches to interpolation
                // (the isosurface sits at d = shell_t, where ramp = 1).
                const double ramp =
                    std::min(1.0, (band - d) / (band - 2.0 * shell_t));
                s = std::min(s, d - shell_t * ramp);
            }
        }
        return (float) s;
    };

    // Two-level fill: exact signed distances on a 4x-coarser lattice, then the
    // fine lattice evaluates exactly only inside a band around the surface and
    // interpolates elsewhere. Distance is 1-Lipschitz, so a fine corner whose
    // interpolated |d| exceeds the band (> half the coarse cell diagonal plus
    // the band the isosurface needs) provably has the interpolation's sign.
    const int C = 4;
    int cdims[3];
    for (int k = 0; k < 3; ++k) cdims[k] = (dims[k] - 2 + C) / C + 1;
    std::vector<float> coarse((size_t) cdims[0] * cdims[1] * cdims[2]);
    parallel_for(coarse.size(), 512, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            const int x = (int)(i % cdims[0]);
            const int y = (int)((i / cdims[0]) % cdims[1]);
            const int z = (int)(i / ((size_t) cdims[0] * cdims[1]));
            coarse[i] = sample({glo[0] + (double) x*C*h, glo[1] + (double) y*C*h, glo[2] + (double) z*C*h});
        }
    });

    // The field is stored negated (inside-high) for mc::extract.
    std::vector<float> field(total);
    parallel_for(total, 4096, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            const int x = (int)(i % dims[0]);
            const int y = (int)((i / dims[0]) % dims[1]);
            const int z = (int)(i / ((size_t) dims[0] * dims[1]));
            const int cx = std::min(x / C, cdims[0] - 2);
            const int cy = std::min(y / C, cdims[1] - 2);
            const int cz = std::min(z / C, cdims[2] - 2);
            const double tx = (double) (x - cx*C) / C, ty = (double) (y - cy*C) / C, tz = (double) (z - cz*C) / C;
            auto cv = [&](int dx, int dy, int dz) {
                return (double) coarse[(size_t)(cx+dx) + (size_t) cdims[0] * ((size_t)(cy+dy) + (size_t) cdims[1] * (cz+dz))];
            };
            const double d0 = (cv(0,0,0)*(1-tx) + cv(1,0,0)*tx) * (1-ty) + (cv(0,1,0)*(1-tx) + cv(1,1,0)*tx) * ty;
            const double d1 = (cv(0,0,1)*(1-tx) + cv(1,0,1)*tx) * (1-ty) + (cv(0,1,1)*(1-tx) + cv(1,1,1)*tx) * ty;
            const double est = d0 * (1-tz) + d1 * tz;
            if (std::fabs(est) >= band) {
                field[i] = (float) -est;
            } else {
                field[i] = -sample({glo[0] + (double) x*h, glo[1] + (double) y*h, glo[2] + (double) z*h});
            }
        }
    });

    // Weld by grid edge: the padded field keeps the surface strictly interior,
    // so exact edge identity dedup applies and near-touching sheets can never
    // fuse into pinholes the way position-quantized welding allows.
    mc::Mesh mesh = mc::extract(field.data(), dims[0], dims[1], dims[2], 0.0f,
                                -1e30f, /*weld_by_grid_edge=*/true);
    field.clear(); field.shrink_to_fit();
    if (mesh.tris.empty()) {
        err = "winding remesh produced an empty mesh";
        return false;
    }
    out_verts.resize(mesh.verts.size());
    for (size_t i = 0; i < mesh.verts.size() / 3; ++i)
        for (int k = 0; k < 3; ++k)
            out_verts[3*i + (size_t) k] = (float)(glo[k] + (double) mesh.verts[3*i + (size_t) k] * h);
    out_normals = std::move(mesh.normals);
    out_tris.assign(mesh.tris.begin(), mesh.tris.end());
    drop_small_islands(out_verts, out_normals, out_tris, 0.0005f);
    if (out_tris.empty()) {
        err = "winding remesh produced only speck components";
        return false;
    }
    // Area-weighted face-based vertex normals: the extractor's field-gradient
    // normals are unreliable where the winding sign or shell term transitions,
    // while the combinatorially oriented faces are trustworthy everywhere.
    std::fill(out_normals.begin(), out_normals.end(), 0.0f);
    for (size_t t = 0; t + 2 < out_tris.size(); t += 3) {
        const int32_t ia = out_tris[t], ib = out_tris[t+1], ic = out_tris[t+2];
        const float * a = out_verts.data() + (size_t) ia * 3;
        const float * b = out_verts.data() + (size_t) ib * 3;
        const float * c = out_verts.data() + (size_t) ic * 3;
        const float ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        const float ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
        const float n[3] = {ab[1]*ac[2] - ab[2]*ac[1],
                            ab[2]*ac[0] - ab[0]*ac[2],
                            ab[0]*ac[1] - ab[1]*ac[0]};
        for (const int32_t i : {ia, ib, ic}) {
            out_normals[(size_t) i*3 + 0] += n[0];
            out_normals[(size_t) i*3 + 1] += n[1];
            out_normals[(size_t) i*3 + 2] += n[2];
        }
    }
    for (size_t i = 0; i < out_normals.size() / 3; ++i) {
        float * n = out_normals.data() + i * 3;
        const float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (len > 1e-20f) { n[0] /= len; n[1] /= len; n[2] /= len; }
    }
    return true;
}

bool winding_project_pbr(const std::vector<float> & source_verts,
                         const std::vector<int32_t> & source_tris,
                         const std::vector<float> & source_pbr,
                         const std::vector<float> & query_points,
                         std::vector<float> & out_pbr,
                         std::string & err) {
    out_pbr.clear();
    const size_t source_nv = source_verts.size() / 3;
    if (source_verts.size() < 9 || source_verts.size() % 3 != 0 ||
        source_tris.size() < 3 || source_tris.size() % 3 != 0) {
        err = "empty PBR projection source";
        return false;
    }
    if (source_pbr.size() != source_nv * 6) {
        err = "PBR projection source has no six-channel material";
        return false;
    }
    if (query_points.size() % 3 != 0) {
        err = "PBR projection query array is not xyz-aligned";
        return false;
    }
    if (query_points.empty()) return true;

    std::vector<int32_t> tris;
    tris.reserve(source_tris.size());
    for (size_t t = 0; t < source_tris.size() / 3; ++t) {
        const int32_t a = source_tris[3*t], b = source_tris[3*t+1], c = source_tris[3*t+2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t) a >= source_nv || (size_t) b >= source_nv || (size_t) c >= source_nv) {
            err = "PBR projection triangle index out of range";
            return false;
        }
        if (a == b || b == c || a == c) continue;
        tris.insert(tris.end(), {a, b, c});
    }
    if (tris.empty()) {
        err = "PBR projection source has no non-degenerate triangles";
        return false;
    }

    BVH bvh;
    bvh.build(source_verts, tris);
    const size_t nq = query_points.size() / 3;
    out_pbr.resize(nq * 6);
    parallel_for(nq, 4096, [&](size_t begin, size_t end) {
        for (size_t qi = begin; qi < end; ++qi) {
            const Vec3 q = {query_points[3*qi], query_points[3*qi+1], query_points[3*qi+2]};
            int32_t tri; double bary[3];
            bvh.closest(q, tri, bary);
            const int32_t * ids = tris.data() + (size_t) tri * 3;
            for (int ch = 0; ch < 6; ++ch)
                out_pbr[6*qi + (size_t) ch] = (float)(
                    bary[0] * source_pbr[6*(size_t) ids[0] + (size_t) ch] +
                    bary[1] * source_pbr[6*(size_t) ids[1] + (size_t) ch] +
                    bary[2] * source_pbr[6*(size_t) ids[2] + (size_t) ch]);
        }
    });
    return true;
}

} // namespace t2print
