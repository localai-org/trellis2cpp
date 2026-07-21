#pragma once

// Print remeshing: turn TRELLIS' arbitrary triangle soup into a closed solid.
// Two backends share this API:
//   - winding_remesh: dependency-free winding-number signed-distance grid +
//     marching tetrahedra. Always compiled; tracks the source surface exactly
//     (no offset) at grid resolution.
//   - alpha_wrap: optional CGAL Alpha Wrap. Builds without CGAL compile the
//     same API and report the backend as unavailable via alpha_available().

#include <cstdint>
#include <string>
#include <vector>

namespace t2print {

// True when any print-remesh backend can run (always, since winding_remesh
// carries no dependency). alpha_available() reports the CGAL backend alone.
bool available();
bool alpha_available();

// alpha_ratio and offset_ratio are fractions of the input bounding-box
// diagonal (for example 0.01 and 0.01/30).  Alpha controls the smallest holes
// and cavities the wrap can enter; offset controls how tightly it encloses the
// source.  Alpha Wrap creates an entirely new surface, so material transfer is
// a separate closest-surface projection step (project_pbr below).
bool alpha_wrap(const std::vector<float> & source_verts,
                const std::vector<int32_t> & source_tris,
                float alpha_ratio,
                float offset_ratio,
                std::vector<float> & out_verts,
                std::vector<float> & out_normals,
                std::vector<int32_t> & out_tris,
                std::string & err);

// Winding-number remesh: sample a signed distance (sign from the generalized
// winding number) on a grid of `grid_res` cells along the longest bounding-box
// axis and contour it.  Closed regions of the soup keep their exact surface;
// open sheets — which enclose no volume — are preserved as thin solid slabs of
// half-thickness max(shell_ratio * diagonal, ~one voxel).  shell_ratio < 0
// disables that fallback (open sheets are then dropped).  Like alpha_wrap the
// output is new geometry, so materials transfer via project_pbr.
bool winding_remesh(const std::vector<float> & source_verts,
                    const std::vector<int32_t> & source_tris,
                    int grid_res,
                    float shell_ratio,
                    std::vector<float> & out_verts,
                    std::vector<float> & out_normals,
                    std::vector<int32_t> & out_tris,
                    std::string & err);

// Project query points onto the closest source triangles and barycentrically
// interpolate their six-channel per-vertex PBR values.  This is the portable
// CPU counterpart of upstream's cuBVH texture reprojection.  Query/output order
// is preserved; query_points contains xyz triples and out_pbr receives six
// floats per query.  Uses CGAL's AABB tree when compiled in, and the
// dependency-free BVH below otherwise.
bool project_pbr(const std::vector<float> & source_verts,
                 const std::vector<int32_t> & source_tris,
                 const std::vector<float> & source_pbr,
                 const std::vector<float> & query_points,
                 std::vector<float> & out_pbr,
                 std::string & err);

// The dependency-free implementation backing project_pbr in non-CGAL builds;
// callable directly for A/B comparison.
bool winding_project_pbr(const std::vector<float> & source_verts,
                         const std::vector<int32_t> & source_tris,
                         const std::vector<float> & source_pbr,
                         const std::vector<float> & query_points,
                         std::vector<float> & out_pbr,
                         std::string & err);

} // namespace t2print
