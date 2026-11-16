#ifndef TTYPT_AREOLE_H
#define TTYPT_AREOLE_H

/*
 * ttypt/areole.h — Generic Distance-Adaptive Multi-Scale Spatial Quadtree & LOD Engine.
 *
 * LICENSE: BSD-2-Clause
 * Copyright (c) 2026, tty-pt
 *
 * Standalone, pure C (C99/C11), domain-agnostic spatial subdivision library.
 * Compatible with native Linux x86_64 (.so) and 32-bit WebAssembly (WASM).
 *
 * Key Architectural Invariants:
 *   1. Configurable Resolution & Pool Arena: Every node spans (base_w * 2^L) x (base_h * 2^L)
 *      world units and samples a configurable (grid_dim + 1) x (grid_dim + 1) lattice allocated
 *      in a single contiguous arena on creation.
 *   2. Bounded Complexity: Quadtree frustum traversal guarantees <= max_active
 *      active leaf nodes across all zoom levels -> O(1) active nodes.
 *   3. Zero Frame Allocations: Fixed pre-allocated arena of node buffers recycled via
 *      LRU tick tracking. Zero malloc/free during real-time camera motion.
 *   4. Pluggable Sampler & Bounds Envelope: Decoupled from any particular noise, game, or
 *      planetary coordinate system.
 *   5. Standard Right-Handed 3D Math: Column-major 4x4 matrices with right-handed LookAt
 *      (s = normalize(f x up), u = s x f) and Gribb-Hartmann 6-plane frustum extraction.
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ttypt/islet_quad.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Canonical Node & Sample Flags */
#define AR_NODE_FLAG_USER0       (1u << 0)
#define AR_NODE_FLAG_DEFERRED    (1u << 30)
#define AR_NODE_FLAG_SPARSE      (1u << 31)
#define AR_NODE_USER_FLAGS_MASK  0x3FFFFFFFu

/* Default Configuration Constants */
#define AR_DEFAULT_VALUE_SCALE          16.0f
#define AR_DEFAULT_MAX_ROOT_CELLS       16u
#define AR_DEFAULT_PROBE_ANCESTOR_CAP   4u
#define AR_DEFAULT_PITCH                0.959931f
#define AR_DEFAULT_DIORAMA_PITCH        AR_DEFAULT_PITCH
#define AR_PITCH_TOP_DOWN               (-1.0f)
#define AR_DEFAULT_TARGET_QUAD_PX       1.0f
#define AR_DEFAULT_ORTHO_QUAD_PX        1.25f
#define AR_DEFAULT_SAG_TANGENT_SLOPE    2.0f
#define AR_DEFAULT_MAX_HYSTERESIS       0.45f
#define AR_DEFAULT_VIEWPORT_W_MULT      8
#define AR_DEFAULT_VIEWPORT_H_MULT      6
#define AR_DEFAULT_CAM_DIST_MULT        2
#define AR_DEFAULT_TOPDOWN_DIST_MULT    4
#define AR_DEFAULT_PERSP_NEAR_RATIO     0.001f
#define AR_DEFAULT_PERSP_FAR_DIST_MULT  10.0f
#define AR_DEFAULT_PERSP_FAR_ZSPAN_MULT 12.0f
#define AR_DEFAULT_MAX_COORD_BOUND      100000000.0f
#define AR_DEFAULT_HORIZON_REACH_MULT   2.5f
#define AR_DEFAULT_MAX_RAY_SPAN_MULT    6.0f
#define AR_DEFAULT_PERSP_Z_LOW_MULT     0.45f
#define AR_DEFAULT_PERSP_Z_HIGH_MULT    0.35f
#define AR_DEFAULT_ORTHO_RAY_T_MIN      (-0.2f)
#define AR_DEFAULT_UNIT_RAY_T_MAX       1.2f
#define AR_DEFAULT_FAR_MID_RAY_T_MAX    1000000.0f
#define AR_DEFAULT_FALLBACK_PERSP_ASPECT 1.5f
#define AR_DEFAULT_BEHIND_CAM_CLIP_W    1.0f
#define AR_DEFAULT_BEHIND_CAM_SCREEN_MULT 2.0f

#define AR_SAMPLE_ATTR_CAP              8u

/* Probe interpolation curve mode */
typedef enum {
	AR_INTERP_SMOOTHSTEP = 0, /* Default C1 cubic Hermite 3t^2 - 2t^3 */
	AR_INTERP_LINEAR     = 1, /* Bilinear u, v */
	AR_INTERP_QUINTIC    = 2, /* C2 quintic 6t^5 - 15t^4 + 10t^3 */
	AR_INTERP_NEAREST    = 3  /* Nearest-neighbor step at 0.5 */
} ar_interp_mode_t;

/* Optional custom 1D interpolation curve callback: maps t in [0, 1] -> weight in [0, 1] */
typedef float (*ar_interp_fn)(float t);

/* Generic sample data evaluated at each grid vertex */
typedef struct {
	int16_t  val;                      /* Generic scalar displacement / field value */
	uint32_t rgba;                     /* Packed 32-bit color or RGBA attributes */
	uint32_t flags;                    /* Generic bitmask flags (bits 0..30 user-defined; bit 31 reserved for AR_NODE_FLAG_SPARSE) */
	uint8_t  attr[AR_SAMPLE_ATTR_CAP]; /* Generic application attribute payload */
} ar_sample_t;

/* User-supplied procedural or discrete spatial sampler */
typedef void (*ar_sampler_fn)(void *user_ctx, uint16_t world, uint8_t level,
                              int32_t wx, int32_t wy, ar_sample_t *out);

/* Optional user-supplied row-batch sampler callback for scanline vectorization */
typedef void (*ar_row_sampler_fn)(void *user_ctx, uint16_t world, uint8_t level,
                                  int32_t wy, int32_t wx_start, int32_t wx_step,
                                  uint32_t count, ar_sample_t *out_row);

/* Optional user-supplied sparsity callback: returns true if discrete data exists */
typedef bool (*ar_sparsity_fn)(void *user_ctx, uint16_t world, uint8_t level,
                               int32_t ox, int32_t oy, int32_t w, int32_t h);

/* Optional user-supplied vertical bounds / envelope callback for curved or displaced surfaces */
typedef void (*ar_bounds_fn)(void *bounds_ctx, float ox, float oy, float w, float h,
                             float *inout_z_min, float *inout_z_max);

/* Optional user-supplied 3D manifold bounding box callback for non-Euclidean / 3D planetary surfaces:
 * Maps 2D manifold patch [ox, ox + w] x [oy, oy + h] with normal displacement [z_min, z_max]
 * into its true 3D Cartesian AABB [out_min[3], out_max[3]], representative 3D center `out_center[3]`,
 * and 3D physical chord span `(out_span_w, out_span_h)` in the coordinate space of `mvp`. */
typedef void (*ar_bounds_3d_fn)(void *bounds_ctx, float ox, float oy, float w, float h,
                                float z_min, float z_max,
                                float out_min[3], float out_max[3],
                                float out_center[3],
                                float *out_span_w, float *out_span_h);

/* Runtime configuration for ar_system_create_ex */
typedef struct {
	int32_t          base_w;               /* Level 0 width in world units (required > 0) */
	int32_t          base_h;               /* Level 0 height in world units (required > 0) */
	uint16_t         grid_dim;             /* Quads per edge (required > 0) */
	uint16_t         pool_capacity;        /* Pre-allocated node pool size (required > 0) */
	uint16_t         max_active;           /* Max active nodes per query (required > 0) */
	uint8_t          max_levels;           /* Max hierarchy levels (required > 0) */
	uint32_t         default_rgba;         /* Fallback RGBA when sample->rgba == 0 */
	float            z_min;                /* Default vertical lower bound when not overridden in ar_cull_opts_t */
	float            z_max;                /* Default vertical upper bound when not overridden in ar_cull_opts_t */
	float            target_quad_px;       /* Default screen-space quad subdivision threshold in pixels */
	float            margin;               /* Default world-space culling margin */
	float            value_scale;          /* Fixed-point scalar divisor (0 = AR_DEFAULT_VALUE_SCALE = 16.0f) */
	uint16_t         max_root_cells;       /* Coarse root grid cell cap before refinement (0 = AR_DEFAULT_MAX_ROOT_CELLS = 16) */
	uint8_t          probe_ancestor_cap;   /* Max ancestor levels searched in ar_system_probe (0 = AR_DEFAULT_PROBE_ANCESTOR_CAP = 4) */
	int32_t          period_x;             /* Optional toroidal world period in world units along X (0 = auto/unbounded) */
	int32_t          period_y;             /* Optional toroidal world period in world units along Y (0 = auto/unbounded) */
	ar_sampler_fn    sampler;              /* Required vertex sampler callback */
	ar_row_sampler_fn row_sampler;         /* Optional fast scanline/row sampler */
	ar_sparsity_fn   sparsity;             /* Optional discrete sparsity callback */
	void            *user_ctx;             /* Opaque context passed to sampler and sparsity */
	uint32_t         max_populates_per_frame; /* Optional frame populate budget (0 = unlimited / use cull_opts) */
	uint32_t         map_capacity;         /* Optional spatial hash map capacity override (0 = pool_capacity * 2) */
	ar_interp_mode_t interp_mode;          /* Probe interpolation curve mode (0 = AR_INTERP_SMOOTHSTEP) */
	ar_interp_fn     interp_fn;            /* Optional custom 1D interpolation weight callback (overrides interp_mode) */
	float            default_pitch;        /* Default pitch in radians when view->pitch == 0 (0 = AR_DEFAULT_PITCH) */
	float            default_cam_dist;     /* Default camera distance when view->eff_cam_dist == 0 (0 = base_w * 2) */
	float            default_topdown_dist; /* Default top-down eye height when view->eff_cam_dist == 0 (0 = base_h * 4) */
	float            default_viewport_w;   /* Default viewport width when view->viewport_w == 0 (0 = base_w * 8) */
	float            default_viewport_h;   /* Default viewport height when view->viewport_h == 0 (0 = base_h * 6) */
	float            sag_radius;           /* Default spherical sag radius when not overridden in view/cull_opts */
	float            sag_tangent_slope;    /* Default linear tangent slope beyond sag_radius (0 = 2.0f) */
	float            split_hysteresis;     /* Default split/merge hysteresis band when not overridden in view/cull_opts */
	ar_bounds_fn     bounds_fn;            /* Default vertical envelope callback when not overridden in cull_opts */
	ar_bounds_3d_fn  bounds_3d_fn;         /* Default 3D manifold AABB callback when not overridden in cull_opts */
	void            *bounds_ctx;           /* Default opaque context passed to bounds_fn / bounds_3d_fn */
	uint8_t          attr_count;           /* Active sample attribute channels (0 = AR_SAMPLE_ATTR_CAP; > AR_SAMPLE_ATTR_CAP = 0) */
	ar_interp_mode_t attr_interp_mode;     /* Attribute probe mode (0 or AR_INTERP_NEAREST = nearest, AR_INTERP_LINEAR/QUINTIC = continuous) */
	bool             interpolate_attrs;    /* When true and attr_interp_mode == 0, interpolate attrs using AR_INTERP_SMOOTHSTEP */
} ar_config_t;

/* Represents an active node in the hierarchical system (public read-only view) */
typedef struct {
	int32_t            ox, oy;          /* Bottom-left origin in world units */
	int32_t            w, h;            /* Physical world extents: (base_w << level), (base_h << level) */
	uint8_t            level;           /* Hierarchy level (0 = base scale) */
	uint16_t           world;           /* World / layer identifier */
	uint32_t           flags;           /* Aggregated node bitmask */
	int16_t            min_val;         /* Minimum scalar value in node */
	int16_t            max_val;         /* Maximum scalar value in node */
	uint32_t           rev;             /* Content revision counter */
	uint16_t           grid_dim;        /* Runtime grid_dim for this node */
	uint16_t           vert_dim;        /* Runtime vert_dim (grid_dim + 1) */
	uint32_t           vert_count;      /* Runtime vert_count (vert_dim * vert_dim) */

	/* Read-only pointers into the single contiguous arena allocated by ar_system_create_ex */
	const int16_t     *values;          /* [vert_count] scalar values for zero-copy GPU upload */
	const uint32_t    *rgba;            /* [vert_count] packed RGBA values */
	const ar_sample_t *samples;         /* [vert_count] full vertex samples */
} ar_node_t;

static inline bool ar_node_is_sparse(const ar_node_t *node)
{
	return node != NULL && ((node->flags & AR_NODE_FLAG_SPARSE) != 0u);
}

/* 6-plane frustum (normalized planes: A*x + B*y + C*z + D >= 0 is inside) */
typedef struct {
	float planes[6][4]; /* Normalized [A, B, C, D] for Left, Right, Bottom, Top, Near, Far */
} ar_frustum_t;

/* Options for matrix frustum culling, vertical envelope, and LOD splitting */
typedef struct {
	float           z_min, z_max;            /* Vertical bounds envelope */
	float           sag_radius;              /* Optional spherical sag radius R (0 = flat or use sys->sag_radius / bounds_fn) */
	float           sag_tangent_slope;       /* Linear tangent slope beyond sag_radius (0 = sys->sag_tangent_slope or 2.0f) */
	float           focal_x, focal_y;        /* Focal center for distance priority and sag computation */
	float           margin;                  /* World-space padding margin added to node bounds */
	float           target_quad_px;          /* Screen-space quad LOD subdivision threshold in pixels */
	float           viewport_w;              /* Viewport width in pixels */
	float           viewport_h;              /* Viewport height in pixels */
	uint8_t         min_level;               /* Minimum LOD level clamp */
	uint8_t         max_level;               /* Maximum LOD level clamp */
	uint8_t         view_id;                 /* View slot identifier (0..31) for multi-view active mask */
	uint16_t        max_populates_per_query; /* Optional per-query new-chunk populate cap (0 = unlimited) */
	uint16_t        max_root_cells;          /* Optional coarse root grid cell cap override (0 = use sys->max_root_cells) */
	float           split_hysteresis;        /* Optional split/merge hysteresis band (e.g. 0.12f) to prevent movement flicker */
	ar_bounds_fn    bounds_fn;               /* Optional generic vertical envelope callback */
	ar_bounds_3d_fn bounds_3d_fn;            /* Optional generic 3D manifold AABB callback */
	void           *bounds_ctx;              /* Opaque context passed to bounds_fn / bounds_3d_fn */
	float        horizon_reach_mult;      /* Optional perspective sag horizon reach multiplier (0 = 2.5f) */
	float        max_ray_span_mult;       /* Optional perspective ray unprojection span multiplier (0 = 6.0f) */
	float        persp_z_low_mult;        /* Optional perspective ray lower Z-plane multiplier (0 = 0.45f) */
	float        persp_z_high_mult;       /* Optional perspective ray upper Z-plane multiplier (0 = 0.35f) */
	float        max_coord_bound;         /* Optional coordinate clamp bound for unprojected rays (0 = 1e8f) */
	float        ortho_ray_t_min;         /* Optional orthographic ray lower t bound (0 = AR_DEFAULT_ORTHO_RAY_T_MIN = -0.2f) */
	float        unit_ray_t_max;          /* Optional unit-cube ray upper t bound (0 = AR_DEFAULT_UNIT_RAY_T_MAX = 1.2f) */
	float        far_mid_ray_t_max;       /* Optional near-to-mid ray upper t bound (0 = AR_DEFAULT_FAR_MID_RAY_T_MAX = 1e6f) */
	float        fallback_persp_aspect;   /* Optional fallback perspective aspect ratio when mvp[0]==0 (0 = 1.5f) */
	float        behind_cam_clip_w;       /* Optional clip-space W threshold for behind/straddling camera (0 = 1.0f) */
	float        behind_cam_screen_mult;  /* Optional viewport width multiplier when straddling camera (0 = 2.0f) */
} ar_cull_opts_t;

/* Camera viewpoint and frustum descriptor */
typedef struct {
	float   cam_x, cam_y, cam_z;     /* Camera focal position in world units */
	float   pitch;                   /* Pitch in radians (0 = default_pitch ~55 deg; AR_PITCH_TOP_DOWN < 0 = top-down) */
	float   yaw;                     /* Yaw in radians */
	float   fov;                     /* Field of view in radians (0 = orthographic) */
	float   zoom;                    /* Global zoom multiplier (>0.0f) */
	float   viewport_w, viewport_h;  /* Viewport dimensions in pixels */
	float   eff_cam_dist;            /* Effective eye-to-target distance */
	uint8_t min_level;               /* Minimum LOD level clamp */
	uint8_t max_level;               /* Maximum LOD level clamp */
	float   near_z;                  /* Optional explicit near clip plane override (0 = auto) */
	float   far_z;                   /* Optional explicit far clip plane override (0 = auto) */
	float   sag_radius;              /* Optional spherical sag radius override (0 = use sys->sag_radius) */
	float   split_hysteresis;        /* Optional split hysteresis override (0 = use sys->split_hysteresis) */
	float   target_quad_px;          /* Optional target quad pixel size override (0 = use sys->target_quad_px) */
	float   margin;                  /* Optional world-space culling margin override (0 = use sys->margin) */
} ar_view_t;

/* Options for 2D orthographic window evaluation (ar_system_eval_ortho_ex) */
typedef struct {
	float   center_x, center_y;      /* World-space center of orthographic window */
	float   world_w, world_h;        /* World-space width and height of orthographic window */
	float   viewport_w, viewport_h;  /* Viewport pixel dimensions (0 = match world_w, world_h) */
	float   target_quad_px;          /* Optional quad pixel threshold override (0 = sys->target_quad_px or 1.25f) */
	uint16_t max_root_cells;         /* Optional max root cells override (0 = sys->max_root_cells or 16) */
	uint8_t min_level;               /* Minimum LOD level clamp */
	uint8_t max_level;               /* Maximum LOD level clamp (0 = sys->max_levels - 1) */
} ar_ortho_opts_t;

/* Multi-View Draw Item returned by ar_system_eval_view / ar_eval_view */
typedef struct {
	const ar_node_t *node;            /* Node data */
	float            world_bounds[4]; /* ox, oy, ox + w, oy + h */
	float            uv_offset[2];    /* UV offset (0.0f, 0.0f for exact match) */
	float            uv_scale;        /* 1.0f for exact match, 0.5f for parent quadrant */
	uint8_t          view_lod;        /* Effective LOD level */
} ar_draw_item_t;

/* Opaque spatial hierarchy engine state */
typedef struct ar_system ar_system_t;
typedef ar_system_t areole_t;

/* Lifecycle Management */
ar_system_t *ar_system_create_ex(const ar_config_t *cfg);

ar_system_t *ar_system_create(int32_t base_w, int32_t base_h,
                              uint16_t grid_dim,
                              uint16_t pool_capacity,
                              uint16_t max_active,
                              uint8_t max_levels,
                              ar_sampler_fn sampler,
                              ar_sparsity_fn sparsity,
                              void *user_ctx);

void ar_system_destroy(ar_system_t *sys);
void ar_system_clear(ar_system_t *sys);
void ar_system_begin_frame(ar_system_t *sys);

size_t ar_system_eval_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                           const ar_view_t *view,
                           ar_draw_item_t *out_items,
                           size_t max_items);

size_t ar_system_eval_ortho(ar_system_t *sys, uint8_t view_id, uint16_t world,
                            float center_x, float center_y,
                            float world_w, float world_h,
                            float viewport_w, float viewport_h);

size_t ar_system_eval_ortho_ex(ar_system_t *sys, uint8_t view_id, uint16_t world,
                               const ar_ortho_opts_t *opts);

float ar_system_get_value_scale(const ar_system_t *sys);
float ar_system_value_to_float(const ar_system_t *sys, int16_t val);
int16_t ar_system_value_from_float(const ar_system_t *sys, float val);
float ar_system_level_texel_span(const ar_system_t *sys, uint8_t level);

/*
 * Matrix Utility Functions (Zero-Dependency Pure C):
 * Column-major 4x4 matrix mathematics matching standard right-handed WebGL / OpenGL conventions.
 *   - ar_mat4_lookat: Standard right-handed LookAt (s = normalize(f x up), u = s x f).
 *   - ar_mat4_lookat_lh: Left-handed / East-flipped LookAt (x = -(up x z), y = x x z).
 */
void ar_mat4_identity(float out[16]);
void ar_mat4_multiply(float out[16], const float a[16], const float b[16]);
void ar_mat4_ortho(float out[16], float left, float right, float bottom, float top, float near_z, float far_z);
void ar_mat4_perspective(float out[16], float fov_rad, float aspect, float near_z, float far_z);
void ar_mat4_lookat(float out[16], const float eye[3], const float target[3], const float up[3]);
void ar_mat4_lookat_lh(float out[16], const float eye[3], const float target[3], const float up[3]);
bool ar_mat4_invert(float out[16], const float m[16]);

/* Frustum Extraction & Testing (Gribb-Hartmann 6-Plane Extraction) */
void ar_frustum_from_matrix(ar_frustum_t *out_frustum, const float mvp[16]);
bool ar_frustum_test_box(const ar_frustum_t *frustum,
                         float min_x, float min_y, float min_z,
                         float max_x, float max_y, float max_z);

/* Primary Matrix Query API */
size_t ar_system_query_matrix(ar_system_t *sys, uint16_t world,
                              const float mvp[16],
                              const ar_cull_opts_t *opts,
                              const ar_node_t **out_nodes,
                              size_t max_out);

size_t ar_system_query_modelview_projection(ar_system_t *sys, uint16_t world,
                                            const float modelview[16],
                                            const float projection[16],
                                            const ar_cull_opts_t *opts,
                                            const ar_node_t **out_nodes,
                                            size_t max_out);

size_t ar_system_update(ar_system_t *sys, uint16_t world,
                        const ar_view_t *view,
                        const ar_node_t **out_nodes,
                        size_t max_out);

bool ar_system_probe(ar_system_t *sys, uint16_t world,
                     int32_t wx, int32_t wy,
                     ar_sample_t *out);

bool ar_system_probe_ex(ar_system_t *sys, uint16_t world,
                        int32_t wx, int32_t wy,
                        uint8_t max_ancestor_levels,
                        ar_sample_t *out);

bool ar_system_probe_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                          int32_t wx, int32_t wy,
                          ar_sample_t *out);

const ar_node_t *ar_system_find(ar_system_t *sys, uint16_t world,
                                uint8_t level, int32_t ox, int32_t oy);

const ar_node_t *ar_system_find_by_key(ar_system_t *sys, uint64_t morton_key);

const ar_node_t *ar_system_ensure(ar_system_t *sys, uint16_t world,
                                  uint8_t level, int32_t ox, int32_t oy);

size_t ar_system_get_active_count(const ar_system_t *sys);
const ar_node_t *ar_system_get_active_node(const ar_system_t *sys, size_t index);

/* Public Node Metadata Accessors */
uint64_t ar_system_node_morton_key(const ar_system_t *sys, const ar_node_t *node);
uint32_t ar_system_node_active_views(const ar_system_t *sys, const ar_node_t *node);
bool ar_system_node_is_active_in_view(const ar_system_t *sys, const ar_node_t *node, uint8_t view_id);

/* Ergonomic Public API Shorthands */
static inline const ar_node_t *ar_node_find(ar_system_t *sys, uint16_t world,
                                            uint8_t level, int32_t ox, int32_t oy)
{
	return ar_system_find(sys, world, level, ox, oy);
}

static inline ar_system_t *ar_create(int32_t base_w, int32_t base_h,
                                     uint16_t grid_dim,
                                     uint16_t pool_capacity,
                                     uint16_t max_active,
                                     uint8_t max_levels,
                                     ar_sampler_fn sampler,
                                     void *user_ctx)
{
	return ar_system_create(base_w, base_h, grid_dim, pool_capacity, max_active, max_levels, sampler, NULL, user_ctx);
}

static inline void ar_destroy(ar_system_t *sys)
{
	ar_system_destroy(sys);
}

static inline void ar_clear(ar_system_t *sys)
{
	ar_system_clear(sys);
}

static inline void ar_begin_frame(ar_system_t *sys)
{
	ar_system_begin_frame(sys);
}

static inline size_t ar_eval_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                                  const ar_view_t *view,
                                  ar_draw_item_t *out_items,
                                  size_t max_items)
{
	return ar_system_eval_view(sys, view_id, world, view, out_items, max_items);
}

static inline size_t ar_eval_ortho(ar_system_t *sys, uint8_t view_id, uint16_t world,
                                   float center_x, float center_y,
                                   float world_w, float world_h,
                                   float viewport_w, float viewport_h)
{
	return ar_system_eval_ortho(sys, view_id, world, center_x, center_y,
	                            world_w, world_h, viewport_w, viewport_h);
}

static inline size_t ar_eval_ortho_ex(ar_system_t *sys, uint8_t view_id, uint16_t world,
                                      const ar_ortho_opts_t *opts)
{
	return ar_system_eval_ortho_ex(sys, view_id, world, opts);
}

static inline size_t ar_update(ar_system_t *sys, uint16_t world,
                               const ar_view_t *view,
                               const ar_node_t **out_nodes,
                               size_t max_out)
{
	return ar_system_update(sys, world, view, out_nodes, max_out);
}

static inline bool ar_probe(ar_system_t *sys, uint16_t world,
                            int32_t wx, int32_t wy,
                            ar_sample_t *out)
{
	return ar_system_probe(sys, world, wx, wy, out);
}

static inline bool ar_probe_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                                 int32_t wx, int32_t wy,
                                 ar_sample_t *out)
{
	return ar_system_probe_view(sys, view_id, world, wx, wy, out);
}

static inline float ar_value_to_float_ex(int16_t val, float scale)
{
	float s = (scale > 0.0f) ? scale : AR_DEFAULT_VALUE_SCALE;
	return (float)val / s;
}

static inline int16_t ar_value_from_float_ex(float val, float scale)
{
	float s = (scale > 0.0f) ? scale : AR_DEFAULT_VALUE_SCALE;
	return (int16_t)(val * s);
}

static inline float ar_value_to_float(int16_t val)
{
	return ar_value_to_float_ex(val, AR_DEFAULT_VALUE_SCALE);
}

static inline int16_t ar_value_from_float(float val)
{
	return ar_value_from_float_ex(val, AR_DEFAULT_VALUE_SCALE);
}

/**
 * Compute the world-space sampling step (texel span) for a quadtree LOD level.
 * Level 0 spans base_span world units across grid_dim cells; each coarser level L
 * doubles the physical span: (base_span << L) / grid_dim.
 */
static inline float ar_level_texel_span(int32_t base_span, uint16_t grid_dim, uint8_t level)
{
	if (grid_dim == 0 || base_span <= 0 || level >= 31)
		return 0.0f;
	float chunk_span = (float)((int64_t)base_span << level);
	return chunk_span / (float)grid_dim;
}

/**
 * Compute the Q16.16 fixed-point world-space sampling step (texel span) for a quadtree LOD level.
 */
static inline int32_t ar_level_texel_span_q16(int32_t base_span, uint16_t grid_dim, uint8_t level)
{
	if (grid_dim == 0 || base_span <= 0 || level >= 31)
		return 0;
	int64_t chunk_span_q16 = ((int64_t)base_span << (level + 16));
	int64_t step_q16 = chunk_span_q16 / (int64_t)grid_dim;
	if (step_q16 > 0x7FFFFFFFLL)
		return 0x7FFFFFFF;
	return (int32_t)step_q16;
}

/* ========================================================================= */
/* Universal Projected Solid-Angle Visibility & Scale-Ordered Spatial Cells  */
/* ========================================================================= */

/**
 * Frame-level camera projection parameters for solid-angle extent culling.
 * Maps world-space bounding extent D = max(w, h) at 3D distance r to projected
 * screen-space pixel size s_px(D, r) = (D * focal_proj) / max(r, min_dist).
 */
typedef struct {
	float focal_proj; /* Pixels per world unit at unit distance (f_proj) */
	float min_px;     /* Minimum projected pixel extent before cutoff (s_min) */
	float full_px;    /* Projected pixel extent for 100% opacity (s_full) */
	float min_dist;   /* Effective minimum camera distance clamp (r_min) */
} ar_solid_angle_t;

static inline void ar_solid_angle_init(ar_solid_angle_t *sa,
                                       float ref_cam_dist,
                                       float viewport_h, float ref_viewport_h,
                                       float eff_cam_dist,
                                       float min_px, float full_px)
{
	if (!sa) return;
	float vp_scale = (viewport_h > 0.0f && ref_viewport_h > 0.0f) ? (viewport_h / ref_viewport_h) : 1.0f;
	sa->focal_proj = (ref_cam_dist > 0.0f ? ref_cam_dist : 220.0f) * vp_scale;
	sa->min_px = (min_px > 0.0f) ? min_px : 1.25f;
	sa->full_px = (full_px > sa->min_px) ? full_px : (sa->min_px + 1.75f);
	sa->min_dist = (eff_cam_dist > 1.0f) ? eff_cam_dist : 1.0f;
}

static inline float ar_solid_angle_dist_2d(const ar_solid_angle_t *sa, float dx, float dy)
{
	float r_min = sa ? sa->min_dist : 1.0f;
	return sqrtf(dx * dx + dy * dy + r_min * r_min);
}

static inline float ar_solid_angle_px(const ar_solid_angle_t *sa, float extent, float dist)
{
	if (!sa || extent <= 0.0f) return 0.0f;
	float r = (dist > sa->min_dist) ? dist : sa->min_dist;
	return (extent * sa->focal_proj) / r;
}

static inline float ar_solid_angle_min_extent(const ar_solid_angle_t *sa, float dist)
{
	if (!sa || sa->focal_proj <= 0.0f) return 0.0f;
	float r = (dist > sa->min_dist) ? dist : sa->min_dist;
	return (sa->min_px * r) / sa->focal_proj;
}

static inline float ar_solid_angle_max_dist(const ar_solid_angle_t *sa, float extent)
{
	if (!sa || sa->min_px <= 0.0f || extent <= 0.0f) return 0.0f;
	return (extent * sa->focal_proj) / sa->min_px;
}

static inline float ar_solid_angle_alpha(const ar_solid_angle_t *sa, float extent, float dist)
{
	float px = ar_solid_angle_px(sa, extent, dist);
	if (!sa || px <= sa->min_px) return 0.0f;
	if (px >= sa->full_px) return 1.0f;
	float t = (px - sa->min_px) / (sa->full_px - sa->min_px);
	return t * t * (3.0f - 2.0f * t);
}

/**
 * Spatial cell summary header for scale-ordered prefix traversal (zero blind replay).
 * Elements in the cell's arena slice [offset .. offset + count - 1] are sorted by
 * quantized extent descending (D_0 >= D_1 >= ... >= D_{count-1}), with North-to-South
 * painter tie-break (y ascending, then x ascending) within each scale class.
 */
typedef struct {
	uint32_t offset;     /* Start index in backing arena */
	uint16_t count;      /* Total elements in cell */
	float    max_extent; /* Largest element extent D_0 in cell */
	float    min_extent; /* Smallest element extent D_{count-1} in cell */
} ar_scale_cell_t;

typedef void (*ar_scale_elem_key_fn)(const void *elem, float *out_extent, float *out_x, float *out_y);

static inline int ar_scale_painter_cmp(float ext_a, float xa, float ya,
                                       float ext_b, float xb, float yb)
{
	int qa = (int)ext_a;
	int qb = (int)ext_b;
	if (qa != qb)
		return (qa > qb) ? -1 : 1; /* Descending scale order */
	if (ya != yb)
		return (ya < yb) ? -1 : 1; /* North-to-South painter order (y ascending) */
	if (xa != xb)
		return (xa < xb) ? -1 : 1;
	return 0;
}

static inline void ar_scale_cell_build(ar_scale_cell_t *cell, uint32_t offset,
                                       void *elems, uint16_t count, size_t elem_size,
                                       ar_scale_elem_key_fn key_fn)
{
	if (!cell) return;
	cell->offset = offset;
	cell->count = count;
	cell->max_extent = 0.0f;
	cell->min_extent = 0.0f;
	if (!elems || count == 0 || elem_size == 0 || !key_fn)
		return;

	uint8_t *base = (uint8_t *)elems;
	uint8_t tmp[256];
	if (count > 1 && elem_size <= sizeof(tmp)) {
		/* Shell sort: O(N log N) in-place with zero heap allocation */
		for (uint16_t gap = count >> 1; gap > 0; gap >>= 1) {
			for (uint16_t i = gap; i < count; i++) {
				memcpy(tmp, base + (size_t)i * elem_size, elem_size);
				float ext_i = 0.0f, xi = 0.0f, yi = 0.0f;
				key_fn(tmp, &ext_i, &xi, &yi);
				uint16_t j = i;
				while (j >= gap) {
					const void *prev = base + (size_t)(j - gap) * elem_size;
					float ext_p = 0.0f, xp = 0.0f, yp = 0.0f;
					key_fn(prev, &ext_p, &xp, &yp);
					if (ar_scale_painter_cmp(ext_p, xp, yp, ext_i, xi, yi) <= 0)
						break;
					memcpy(base + (size_t)j * elem_size, prev, elem_size);
					j -= gap;
				}
				memcpy(base + (size_t)j * elem_size, tmp, elem_size);
			}
		}
	}

	float max_e = 0.0f, min_e = 1e30f;
	for (uint16_t i = 0; i < count; i++) {
		float e = 0.0f, x = 0.0f, y = 0.0f;
		key_fn(base + (size_t)i * elem_size, &e, &x, &y);
		if (e > max_e) max_e = e;
		if (e < min_e) min_e = e;
	}
	cell->max_extent = max_e;
	cell->min_extent = (min_e <= max_e) ? min_e : 0.0f;
}

static inline bool ar_scale_cell_cull(const ar_scale_cell_t *cell,
                                      const ar_solid_angle_t *sa,
                                      float cell_dist, float cell_radius,
                                      float *out_cutoff_extent)
{
	if (!cell || cell->count == 0 || !sa)
		return false;
	float r_near = cell_dist - cell_radius;
	if (r_near < sa->min_dist)
		r_near = sa->min_dist;
	float near_cutoff = ar_solid_angle_min_extent(sa, r_near);
	if (cell->max_extent < near_cutoff)
		return false;
	if (out_cutoff_extent)
		*out_cutoff_extent = ar_solid_angle_min_extent(sa, cell_dist);
	return true;
}

/* Statistics inspection for diagnostics and unit tests */
typedef struct {
	size_t  pool_capacity;
	size_t  pool_in_use;
	size_t  active_chunks_last_frame;
	size_t  sampler_invocations_total;
	size_t  cache_hits_total;
	size_t  cache_misses_total;
	uint8_t max_level_observed;
	uint8_t min_level_observed;
} ar_stats_t;

void ar_system_get_stats(const ar_system_t *sys, ar_stats_t *out_stats);

#ifdef __cplusplus
}
#endif

#endif /* TTYPT_AREOLE_H */
