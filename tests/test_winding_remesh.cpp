// Winding-number print remesh (winding_remesh.cpp). Unlike the CGAL test this
// backend is always compiled, so the suite never skips.

#include "mesh_export.h"
#include "print_remesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <vector>

struct EdgeUse { int count = 0; int direction = 0; };

// Closed + 2-manifold + consistently wound, and the signed volume (divergence
// theorem; positive iff wound outward).
static bool check_solid(const std::vector<float> & verts,
                        const std::vector<int32_t> & tris,
                        const char * label, double & volume) {
    std::unordered_map<uint64_t, EdgeUse> edges;
    double signed_volume_6 = 0.0;
    for (size_t t = 0; t < tris.size(); t += 3) {
        const int32_t a = tris[t], b = tris[t+1], c = tris[t+2];
        for (const auto e : {std::pair<int32_t,int32_t>{a,b}, {b,c}, {c,a}}) {
            const uint32_t lo = (uint32_t) std::min(e.first, e.second);
            const uint32_t hi = (uint32_t) std::max(e.first, e.second);
            EdgeUse & use = edges[((uint64_t) lo << 32) | hi];
            use.count++;
            use.direction += e.first < e.second ? 1 : -1;
        }
        const float * va = verts.data() + (size_t) a * 3;
        const float * vb = verts.data() + (size_t) b * 3;
        const float * vc = verts.data() + (size_t) c * 3;
        signed_volume_6 += va[0] * (vb[1]*vc[2] - vb[2]*vc[1])
                         + va[1] * (vb[2]*vc[0] - vb[0]*vc[2])
                         + va[2] * (vb[0]*vc[1] - vb[1]*vc[0]);
    }
    for (const auto & it : edges) {
        if (it.second.count != 2 || it.second.direction != 0) {
            std::fprintf(stderr, "%s: non-manifold/inconsistently wound edge: count=%d direction=%d\n",
                         label, it.second.count, it.second.direction);
            return false;
        }
    }
    volume = signed_volume_6 / 6.0;
    if (volume <= 0.0) {
        std::fprintf(stderr, "%s: enclosed volume %.6g is not positive (inward winding?)\n",
                     label, volume);
        return false;
    }
    return true;
}

static void append_box(std::vector<float> & verts, std::vector<int32_t> & tris,
                       float x0, float y0, float z0, float x1, float y1, float z1) {
    const int32_t base = (int32_t)(verts.size() / 3);
    const float corner[8][3] = {
        {x0,y0,z0},{x1,y0,z0},{x1,y1,z0},{x0,y1,z0},
        {x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1},
    };
    for (const auto & c : corner) verts.insert(verts.end(), {c[0], c[1], c[2]});
    // Outward-wound faces (CCW seen from outside).
    const int32_t face[12][3] = {
        {0,2,1},{0,3,2},           // z = z0 (normal -z)
        {4,5,6},{4,6,7},           // z = z1 (normal +z)
        {0,1,5},{0,5,4},           // y = y0 (normal -y)
        {3,6,2},{3,7,6},           // y = y1 (normal +y)
        {0,4,7},{0,7,3},           // x = x0 (normal -x)
        {1,2,6},{1,6,5},           // x = x1 (normal +x)
    };
    for (const auto & f : face) tris.insert(tris.end(), {base+f[0], base+f[1], base+f[2]});
}

int main() {
    std::string err;

    // 1) An open zero-thickness triangle is exactly the winding number's blind
    //    spot (it encloses no volume). The shell fallback must still produce a
    //    closed, printable solid around it.
    {
        const std::vector<float> verts = {0,0,0, 1,0,0, 0,1,0};
        const std::vector<int32_t> tris = {0,1,2};
        std::vector<float> ov, on;
        std::vector<int32_t> ot;
        if (!t2print::winding_remesh(verts, tris, 64, 0.0f, ov, on, ot, err)) {
            std::fprintf(stderr, "open-sheet winding_remesh failed: %s\n", err.c_str());
            return 1;
        }
        if (ov.empty() || ot.empty() || on.size() != ov.size()) {
            std::fprintf(stderr, "open-sheet output incomplete: %zu verts %zu tris %zu normals\n",
                         ov.size()/3, ot.size()/3, on.size()/3);
            return 1;
        }
        double volume = 0.0;
        if (!check_solid(ov, ot, "open-sheet", volume)) return 1;
    }

    // 2) Two overlapping outward-wound cubes: a self-intersecting soup whose
    //    union is well defined. The remesh must reproduce that union without
    //    the shell dilating the closed surfaces: vol = 2 - 0.5^3 = 1.875.
    {
        std::vector<float> verts;
        std::vector<int32_t> tris;
        append_box(verts, tris, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);
        append_box(verts, tris, 0.5f, 0.5f, 0.5f, 1.5f, 1.5f, 1.5f);
        std::vector<float> ov, on;
        std::vector<int32_t> ot;
        if (!t2print::winding_remesh(verts, tris, 96, 0.0f, ov, on, ot, err)) {
            std::fprintf(stderr, "cube-union winding_remesh failed: %s\n", err.c_str());
            return 1;
        }
        double volume = 0.0;
        if (!check_solid(ov, ot, "cube-union", volume)) return 1;
        const double expected = 2.0 - 0.125;
        if (std::fabs(volume - expected) > 0.06 * expected) {
            std::fprintf(stderr, "cube-union volume %.4f, expected ~%.4f\n", volume, expected);
            return 1;
        }
        // The surface must track the input: no vertex may sit further than two
        // voxels (~1.5/96 each) outside the exact union's bounding box.
        const double tol = 2.0 * 1.5 / 96.0;
        for (size_t i = 0; i < ov.size() / 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                const double v = ov[3*i + (size_t) k];
                if (v < 0.0 - tol || v > 1.5 + tol) {
                    std::fprintf(stderr, "cube-union vertex escapes the union: axis %d value %.4f\n", k, v);
                    return 1;
                }
            }
        }
    }

    // 3) The full prepare_print_mesh path with the env-var backend override the
    //    A/B prototype exposes (works in CGAL and non-CGAL builds alike).
    {
        setenv("T2_PRINT_BACKEND", "winding", 1);
        setenv("T2_PRINT_RES", "64", 1);
        const float verts[] = {0,0,0, 1,0,0, 0,1,0};
        const int32_t tris[] = {0,1,2};
        const float pbr[] = {
            1,0,0,0,0.5f,1, 0,1,0,0,0.5f,1, 0,0,1,0,0.5f,1,
        };
        t2glb::MeshExportOptions opt;
        opt.components = t2glb::ComponentFilter::KeepAll;
        t2glb::PreparedMesh out;
        if (!t2glb::prepare_print_mesh(verts, 3, tris, 1, pbr, opt,
                                       0.20f, 0.03f, out, err)) {
            std::fprintf(stderr, "prepare_print_mesh(winding) failed: %s\n", err.c_str());
            return 1;
        }
        unsetenv("T2_PRINT_BACKEND");
        unsetenv("T2_PRINT_RES");
        if (out.verts.empty() || out.tris.empty() || out.normals.size() != out.verts.size()) {
            std::fprintf(stderr, "prepare_print_mesh(winding) output incomplete\n");
            return 1;
        }
        // Textured sources get a per-vertex preview material via projection.
        if (out.pbr.size() != out.verts.size() / 3 * 6) {
            std::fprintf(stderr, "prepare_print_mesh(winding) missing projected preview PBR\n");
            return 1;
        }
    }

    // 4) The dependency-free closest-point projection must interpolate
    //    barycentrically, exactly like the CGAL backend it substitutes.
    {
        const std::vector<float> verts = {0,0,0, 1,0,0, 0,1,0};
        const std::vector<int32_t> tris = {0,1,2};
        const std::vector<float> pbr = {
            1,0,0,0,0.5f,1, 0,1,0,0,0.5f,1, 0,0,1,0,0.5f,1,
        };
        const std::vector<float> queries = {0.25f, 0.25f, 0.5f};
        std::vector<float> projected;
        if (!t2print::winding_project_pbr(verts, tris, pbr, queries, projected, err) ||
            projected.size() != 6) {
            std::fprintf(stderr, "winding_project_pbr failed: %s\n", err.c_str());
            return 1;
        }
        const float expected[] = {0.5f, 0.25f, 0.25f, 0.0f, 0.5f, 1.0f};
        for (int i = 0; i < 6; ++i) {
            if (std::fabs(projected[i] - expected[i]) > 1e-5f) {
                std::fprintf(stderr, "winding_project_pbr channel %d: got %.7g expected %.7g\n",
                             i, projected[i], expected[i]);
                return 1;
            }
        }
    }

    std::printf("RESULT: PASS\n");
    return 0;
}
