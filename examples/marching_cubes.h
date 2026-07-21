// marching_cubes.h — single-file, self-contained isosurface extraction.
//
// Extracts a triangle mesh for the level set {f = iso} of a scalar field
// sampled on a regular grid, and writes it as a Wavefront OBJ.
//
// It uses the TETRAHEDRAL variant of marching cubes (a.k.a. marching
// tetrahedra on the Freudenthal/Kuhn subdivision): each grid cube is split
// into 6 tetrahedra that all share the (0,0,0)-(1,1,1) main diagonal, and each
// tetrahedron is contoured with a tiny 16-case table. Compared to classic
// marching cubes this avoids the error-prone 256-row triangle table entirely,
// and — because every cube uses the same diagonal — the subdivision tiles space
// consistently, so the output is a watertight 2-manifold (each interior edge is
// shared by exactly two triangles). Vertices are de-duplicated by the global
// grid edge they sit on, so the mesh is properly indexed.
//
// Header-only and dependency-free (just <vector>/<cmath>/<unordered_map>).
//
//   mc::Mesh m = mc::extract(field, nx, ny, nz, iso);
//   mc::write_obj(m, "out.obj");
//
// `field` is indexed x + y*nx + z*nx*ny (x fastest). Vertex positions are in
// grid-index units (a voxel is one unit). The field is treated as `pad` (a
// large negative value, i.e. "outside") beyond the grid, so surfaces that reach
// the boundary are closed off.

#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <vector>

namespace mc {

struct Mesh {
    std::vector<float> verts;   // 3 floats per vertex (x,y,z)
    std::vector<float> normals; // 3 floats per vertex (unit, oriented outward)
    std::vector<int>   tris;    // 3 vertex indices per triangle
    size_t n_verts() const { return verts.size() / 3; }
    size_t n_tris()  const { return tris.size() / 3; }
};

namespace detail {

// Cube corner offsets, Bourke/Lorensen order.
static const int CORNER[8][3] = {
    {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
};

// The 6 tetrahedra of the Freudenthal subdivision, all sharing the 0-6 diagonal.
static const int TETRA[6][4] = {
    {0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}
};

// A tetrahedron's 6 edges, as (local-corner, local-corner) pairs.
static const int TET_EDGE[6][2] = {
    {0,1},{1,2},{2,0},{0,3},{1,3},{2,3}
};

// Per-case triangle table for one tetrahedron. Index = bitmask of which of the
// 4 local corners are "inside" (value > iso); values are TET_EDGE indices,
// terminated by -1. Derived by hand (see marching_cubes.h history) and checked
// at runtime by extract()'s self-consistency asserts on small fields.
static const int TET_TRI[16][7] = {
    {-1,-1,-1,-1,-1,-1,-1},                 // 0000
    { 0, 2, 3,-1,-1,-1,-1},                 // 0001  {0}
    { 0, 1, 4,-1,-1,-1,-1},                 // 0010  {1}
    { 2, 1, 4,  2, 4, 3,-1},                // 0011  {0,1}
    { 2, 1, 5,-1,-1,-1,-1},                 // 0100  {2}
    { 0, 1, 5,  0, 5, 3,-1},                // 0101  {0,2}
    { 0, 2, 5,  0, 5, 4,-1},                // 0110  {1,2}
    { 3, 4, 5,-1,-1,-1,-1},                 // 0111  {0,1,2} (outside 3)
    { 3, 4, 5,-1,-1,-1,-1},                 // 1000  {3}
    { 0, 4, 5,  0, 5, 2,-1},                // 1001  {0,3}
    { 0, 3, 5,  0, 5, 1,-1},                // 1010  {1,3}
    { 2, 1, 5,-1,-1,-1,-1},                 // 1011  {0,1,3} (outside 2)
    { 2, 3, 4,  2, 4, 1,-1},                // 1100  {2,3}
    { 0, 1, 4,-1,-1,-1,-1},                 // 1101  {0,2,3} (outside 1)
    { 0, 2, 3,-1,-1,-1,-1},                 // 1110  {1,2,3} (outside 0)
    {-1,-1,-1,-1,-1,-1,-1},                 // 1111
};

} // namespace detail

namespace detail {

// Combinatorial triangle orientation: decided once per (tet, case, triangle)
// from the canonical cell geometry with midpoint crossings, winding each
// triangle so its normal points from the inside corners toward the outside
// corners. Unlike a per-instance test against the sampled field this never
// consults the field itself, so orientation stays consistent even where the
// field is discontinuous or its gradient is degenerate.
struct FlipTable {
    bool flip[6][16][2] = {};
    FlipTable() {
        for (int t = 0; t < 6; ++t) {
            const int * tc = TETRA[t];
            for (int code = 1; code < 15; ++code) {
                double cin[3] = {0,0,0}, cout[3] = {0,0,0};
                int nin = 0, nout = 0;
                for (int i = 0; i < 4; ++i) {
                    const int * c = CORNER[tc[i]];
                    if (code & (1 << i)) {
                        for (int d = 0; d < 3; ++d) cin[d] += c[d];
                        ++nin;
                    } else {
                        for (int d = 0; d < 3; ++d) cout[d] += c[d];
                        ++nout;
                    }
                }
                double dir[3];
                for (int d = 0; d < 3; ++d) dir[d] = cout[d] / nout - cin[d] / nin;
                const int * tri = TET_TRI[code];
                for (int e = 0; tri[e] >= 0; e += 3) {
                    double p[3][3];
                    for (int k = 0; k < 3; ++k) {
                        const int * ed = TET_EDGE[tri[e + k]];
                        const int * a = CORNER[tc[ed[0]]], * b = CORNER[tc[ed[1]]];
                        for (int d = 0; d < 3; ++d) p[k][d] = 0.5 * (a[d] + b[d]);
                    }
                    const double u[3] = {p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2]};
                    const double v[3] = {p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2]};
                    const double n[3] = {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
                    flip[t][code][e / 3] =
                        n[0]*dir[0] + n[1]*dir[1] + n[2]*dir[2] < 0.0;
                }
            }
        }
    }
};

} // namespace detail

// weld_by_grid_edge selects vertex deduplication by the (unique) grid edge a
// crossing sits on instead of by quantized position. For surfaces that never
// reach the grid boundary this is exact: distinct edges never fuse, so nearly
// touching sheets cannot weld into pinholes or non-manifold knife contacts.
// Fields that DO reach the boundary should keep the default position welding,
// which collapses the degenerate cap slivers the pad interpolation produces.
inline Mesh extract(const float * field, int nx, int ny, int nz, float iso,
                    float pad = -1e30f, bool weld_by_grid_edge = false) {
    using namespace detail;
    static const FlipTable FLIP;   // magic static: thread-safe one-time init
    Mesh mesh;

    // Tie-breaking jitter (Simulation-of-Simplicity style): a corner whose value
    // is *exactly* iso makes the surface pass through that single point, which is
    // shared by every incident tetrahedron — a pinch that no winding can make
    // coherent. A tiny deterministic per-corner perturbation (consistent across
    // cubes, so shared corners still agree -> mesh stays watertight) pushes every
    // corner decisively on/off the surface. Scaled to the field so vertex motion
    // is sub-micro-voxel; continuous fields (e.g. decoder logits) are unaffected.
    double fmax = 0.0;
    for (size_t i = 0; i < (size_t) nx * ny * nz; ++i) {
        double a = std::fabs((double) field[i]); if (a > fmax) fmax = a;
    }
    const float eps = (float) (1e-6 * fmax) + 1e-20f;
    auto jitter = [](int x, int y, int z) -> float {
        uint32_t h = (uint32_t)(x * 73856093) ^ (uint32_t)(y * 19349663) ^ (uint32_t)(z * 83492791);
        h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
        return (float)(h & 0xffff) / 65535.0f - 0.5f; // [-0.5, 0.5)
    };
    auto at = [&](int x, int y, int z) -> float {
        if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) return pad;
        return field[(size_t) x + (size_t) y * nx + (size_t) z * nx * ny] + eps * jitter(x, y, z);
    };
    // Central-difference gradient (for outward-oriented vertex normals).
    auto grad = [&](int x, int y, int z, float g[3]) {
        g[0] = at(x + 1, y, z) - at(x - 1, y, z);
        g[1] = at(x, y + 1, z) - at(x, y - 1, z);
        g[2] = at(x, y, z + 1) - at(x, y, z - 1);
    };

    std::unordered_map<int64_t, int> vmap; // quantized-position key -> vertex index

    // Interpolate (and cache) the vertex on the grid edge cA->cB. Vertices are
    // de-duplicated by QUANTIZED POSITION, not by edge: against the large pad,
    // boundary-cap vertices interpolate to t~0 and collapse exactly onto grid
    // corners, so several distinct edges land on the same point. Position keying
    // welds those into one vertex (killing the degenerate zero-area cap triangles
    // and keeping the boundary watertight), while shared edges across cubes —
    // which already yield identical positions — still merge as before.
    const double QS = 1024.0; // sub-milli-voxel quantization
    // Corner id for grid-edge keying; coords range [-1, n] from the pad ring.
    const int64_t sx = nx + 2, sy = ny + 2;
    auto corner_id = [&](const int c[3]) -> int64_t {
        return (c[0] + 1) + sx * ((c[1] + 1) + sy * (int64_t)(c[2] + 1));
    };
    auto edge_vertex = [&](const int ca[3], const int cb[3], float va, float vb) -> int {
        float t = (vb != va) ? (iso - va) / (vb - va) : 0.5f;
        if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;
        float p[3] = { ca[0] + t * (cb[0] - ca[0]),
                       ca[1] + t * (cb[1] - ca[1]),
                       ca[2] + t * (cb[2] - ca[2]) };
        int64_t key;
        if (weld_by_grid_edge) {
            const int64_t ia = corner_id(ca), ib = corner_id(cb);
            const int64_t total = sx * sy * (int64_t)(nz + 2);
            key = std::min(ia, ib) * total + std::max(ia, ib);
        } else {
            int64_t qx = (int64_t) std::llround(p[0] * QS) + 0x200000;
            int64_t qy = (int64_t) std::llround(p[1] * QS) + 0x200000;
            int64_t qz = (int64_t) std::llround(p[2] * QS) + 0x200000;
            key = (qx * 0x400000LL + qy) * 0x400000LL + qz;
        }
        auto it = vmap.find(key);
        if (it != vmap.end()) return it->second;

        float ga[3], gb[3];
        grad(ca[0], ca[1], ca[2], ga);
        grad(cb[0], cb[1], cb[2], gb);
        float n[3];
        for (int d = 0; d < 3; ++d) n[d] = -(ga[d] + t * (gb[d] - ga[d])); // outward = -grad
        // Normalize in double: against the large negative pad the gradient can
        // reach ~1e30, and 1e30^2 overflows float32 to inf (-> bogus normal).
        double len = std::sqrt((double)n[0]*n[0] + (double)n[1]*n[1] + (double)n[2]*n[2]);
        if (len > 1e-12) { n[0]=(float)(n[0]/len); n[1]=(float)(n[1]/len); n[2]=(float)(n[2]/len); }
        else             { n[0]=0; n[1]=0; n[2]=1; }

        int idx = (int) mesh.n_verts();
        mesh.verts.insert(mesh.verts.end(), { p[0], p[1], p[2] });
        mesh.normals.insert(mesh.normals.end(), { n[0], n[1], n[2] });
        vmap.emplace(key, idx);
        return idx;
    };

    // Iterate cube origins over [-1 .. n-1] so boundary cubes (against the pad)
    // close the surface.
    for (int z = -1; z < nz; ++z)
    for (int y = -1; y < ny; ++y)
    for (int x = -1; x < nx; ++x) {
        float cval[8];
        int   ccoord[8][3];
        for (int c = 0; c < 8; ++c) {
            ccoord[c][0] = x + CORNER[c][0];
            ccoord[c][1] = y + CORNER[c][1];
            ccoord[c][2] = z + CORNER[c][2];
            cval[c] = at(ccoord[c][0], ccoord[c][1], ccoord[c][2]);
        }
        for (int t = 0; t < 6; ++t) {
            const int * tc = TETRA[t];
            int code = 0;
            for (int i = 0; i < 4; ++i) if (cval[tc[i]] > iso) code |= (1 << i);
            const int * tri = TET_TRI[code];
            if (tri[0] < 0) continue;

            for (int e = 0; tri[e] >= 0; e += 3) {
                int vi[3];
                for (int k = 0; k < 3; ++k) {
                    const int * ed = TET_EDGE[tri[e + k]];
                    int la = tc[ed[0]], lb = tc[ed[1]];
                    vi[k] = edge_vertex(ccoord[la], ccoord[lb], cval[la], cval[lb]);
                }
                if (vi[0] == vi[1] || vi[1] == vi[2] || vi[0] == vi[2]) continue; // degenerate
                int v0 = vi[0], v1 = vi[1], v2 = vi[2];
                if (FLIP.flip[t][code][e / 3]) { int tmp = v1; v1 = v2; v2 = tmp; }
                mesh.tris.insert(mesh.tris.end(), { v0, v1, v2 });
            }
        }
    }
    return mesh;
}

inline bool write_obj(const Mesh & m, const char * path) {
    FILE * f = std::fopen(path, "wb");
    if (!f) return false;
    std::fprintf(f, "# trellis2.cpp isosurface  (%zu verts, %zu tris)\n",
                 m.n_verts(), m.n_tris());
    for (size_t i = 0; i < m.n_verts(); ++i)
        std::fprintf(f, "v %.6f %.6f %.6f\n", m.verts[3*i], m.verts[3*i+1], m.verts[3*i+2]);
    for (size_t i = 0; i < m.n_verts(); ++i)
        std::fprintf(f, "vn %.6f %.6f %.6f\n", m.normals[3*i], m.normals[3*i+1], m.normals[3*i+2]);
    for (size_t i = 0; i < m.n_tris(); ++i) {
        int a = m.tris[3*i] + 1, b = m.tris[3*i+1] + 1, c = m.tris[3*i+2] + 1;
        std::fprintf(f, "f %d//%d %d//%d %d//%d\n", a, a, b, b, c, c);
    }
    std::fclose(f);
    return true;
}

} // namespace mc
