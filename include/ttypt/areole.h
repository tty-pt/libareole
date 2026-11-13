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

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MORTON_QUAD_DEFINED
#define MORTON_QUAD_DEFINED

static inline int32_t morton_fdiv(int32_t a, int32_t b)
{
	int32_t q = a / b;
	int32_t r = a % b;
	if (r != 0 && ((a < 0) ^ (b < 0)))
		q--;
	return q;
}

static inline uint32_t morton_spread_dense16(uint16_t v)
{
#if defined(__BMI2__)
	return (uint32_t)__builtin_ia32_pdep_di(v, 0x55555555ULL);
#else
	uint32_t x = v;
	x = (x | (x << 8)) & 0x00FF00FFu;
	x = (x | (x << 4)) & 0x0F0F0F0Fu;
	x = (x | (x << 2)) & 0x33333333u;
	x = (x | (x << 1)) & 0x55555555u;
	return x;
#endif
}

static inline uint16_t morton_compact_dense16(uint32_t v)
{
#if defined(__BMI2__)
	return (uint16_t)__builtin_ia32_pext_di(v, 0x55555555ULL);
#else
	uint32_t x = v & 0x55555555u;
	x = (x | (x >> 1)) & 0x33333333u;
	x = (x | (x >> 2)) & 0x0F0F0F0Fu;
	x = (x | (x >> 4)) & 0x00FF00FFu;
	x = (x | (x >> 8)) & 0x0000FFFFu;
	return (uint16_t)x;
#endif
}

#define MORTON_QUAD_BIAS     32768

static inline uint64_t morton_quad_pack(uint16_t world, uint8_t level, int32_t cell_x, int32_t cell_y)
{
	int32_t clamped_x = cell_x < -32768 ? -32768 : (cell_x > 32767 ? 32767 : cell_x);
	int32_t clamped_y = cell_y < -32768 ? -32768 : (cell_y > 32767 ? 32767 : cell_y);
	uint16_t ux = (uint16_t)(clamped_x + MORTON_QUAD_BIAS);
	uint16_t uy = (uint16_t)(clamped_y + MORTON_QUAD_BIAS);
	uint32_t m = morton_spread_dense16(ux) | (morton_spread_dense16(uy) << 1);
	uint64_t key = ((uint64_t)(world & 0x3FFF) << 40) | ((uint64_t)(level & 0x1F) << 32) | (uint64_t)m;
	return key;
}

static inline void morton_quad_unpack(uint64_t key, uint16_t *world, uint8_t *level, int32_t *cell_x, int32_t *cell_y)
{
	uint8_t lvl = (uint8_t)((key >> 32) & 0x1F);
	uint32_t m = (uint32_t)(key & 0xFFFFFFFFu);
	uint16_t ux = morton_compact_dense16(m);
	uint16_t uy = morton_compact_dense16(m >> 1);
	if (world) *world = (uint16_t)((key >> 40) & 0x3FFF);
	if (level) *level = lvl;
	if (cell_x) *cell_x = (int32_t)ux - MORTON_QUAD_BIAS;
	if (cell_y) *cell_y = (int32_t)uy - MORTON_QUAD_BIAS;
}

static inline uint64_t morton_quad_key(uint16_t world, uint8_t level, int32_t ox, int32_t oy, int32_t base_w, int32_t base_h)
{
	int32_t w = base_w << (level < 16 ? level : 16);
	int32_t h = base_h << (level < 16 ? level : 16);
	int32_t cx = (w > 0) ? morton_fdiv(ox, w) : 0;
	int32_t cy = (h > 0) ? morton_fdiv(oy, h) : 0;
	return morton_quad_pack(world, level, cx, cy);
}

static inline void morton_quad_decode(uint64_t key, int32_t base_w, int32_t base_h,
                                      uint16_t *world, uint8_t *level, int32_t *ox, int32_t *oy)
{
	int32_t cx = 0, cy = 0;
	uint8_t lvl = 0;
	morton_quad_unpack(key, world, &lvl, &cx, &cy);
	if (level) *level = lvl;
	int32_t w = base_w << (lvl < 16 ? lvl : 16);
	int32_t h = base_h << (lvl < 16 ? lvl : 16);
	if (ox) *ox = cx * w;
	if (oy) *oy = cy * h;
}

static inline uint64_t morton_quad_parent(uint64_t key)
{
	uint16_t world = 0;
	uint8_t lvl = 0;
	int32_t cx = 0, cy = 0;
	morton_quad_unpack(key, &world, &lvl, &cx, &cy);
	if (lvl >= 16) return key;
	return morton_quad_pack(world, (uint8_t)(lvl + 1), morton_fdiv(cx, 2), morton_fdiv(cy, 2));
}

static inline uint64_t morton_quad_child(uint64_t key, uint8_t quadrant)
{
	uint16_t world = 0;
	uint8_t lvl = 0;
	int32_t cx = 0, cy = 0;
	morton_quad_unpack(key, &world, &lvl, &cx, &cy);
	if (lvl == 0) return key;
	int32_t child_cx = (cx << 1) | (quadrant & 1);
	int32_t child_cy = (cy << 1) | ((quadrant >> 1) & 1);
	return morton_quad_pack(world, (uint8_t)(lvl - 1), child_cx, child_cy);
}

static inline uint8_t morton_quad_quadrant(uint64_t key)
{
	int32_t cx = 0, cy = 0;
	morton_quad_unpack(key, NULL, NULL, &cx, &cy);
	return (uint8_t)((cx & 1) | ((cy & 1) << 1));
}

#endif /* MORTON_QUAD_DEFINED */

/* Generic sample data evaluated at each grid vertex */
typedef struct {
	int16_t  val;        /* Generic scalar displacement / field value */
	uint32_t rgba;       /* Packed 32-bit color or RGBA attributes */
	uint32_t flags;      /* Generic bitmask flags (e.g. discrete, surface type, etc.) */
	uint8_t  attr[8];    /* Generic 8-byte application payload */
} ar_sample_t;

/* User-supplied procedural or discrete spatial sampler */
typedef void (*ar_sampler_fn)(void *user_ctx, uint16_t world, uint8_t level,
                              int32_t wx, int32_t wy, ar_sample_t *out);

/* Optional user-supplied sparsity callback: returns true if discrete data exists */
typedef bool (*ar_sparsity_fn)(void *user_ctx, uint16_t world, uint8_t level,
                               int32_t ox, int32_t oy, int32_t w, int32_t h);

/* Optional user-supplied vertical bounds / envelope callback for curved or displaced surfaces */
typedef void (*ar_bounds_fn)(void *bounds_ctx, float ox, float oy, float w, float h,
                             float *inout_z_min, float *inout_z_max);

/* Runtime configuration for ar_system_create_ex */
typedef struct {
	int32_t        base_w;          /* Level 0 width in world units (required > 0) */
	int32_t        base_h;          /* Level 0 height in world units (required > 0) */
	uint16_t       grid_dim;        /* Quads per edge (required > 0) */
	uint16_t       pool_capacity;   /* Pre-allocated node pool size (required > 0) */
	uint16_t       max_active;      /* Max active nodes per query (required > 0) */
	uint8_t        max_levels;      /* Max hierarchy levels (required > 0) */
	uint32_t       default_rgba;    /* Fallback RGBA when sample->rgba == 0 */
	float          z_min;           /* Default vertical lower bound when not overridden in ar_cull_opts_t */
	float          z_max;           /* Default vertical upper bound when not overridden in ar_cull_opts_t */
	float          target_quad_px;  /* Default screen-space quad subdivision threshold in pixels */
	float          margin;          /* Default world-space culling margin */
	ar_sampler_fn  sampler;         /* Required vertex sampler callback */
	ar_sparsity_fn sparsity;        /* Optional discrete sparsity callback */
	void          *user_ctx;        /* Opaque context passed to sampler and sparsity */
} ar_config_t;

/* Represents an active node in the hierarchical system */
typedef struct {
	int32_t        ox, oy;          /* Bottom-left origin in world units */
	int32_t        w, h;            /* Physical world extents: (base_w << level), (base_h << level) */
	uint8_t        level;           /* Hierarchy level (0 = base scale) */
	uint16_t       world;           /* World / layer identifier */
	uint32_t       last_used_tick;  /* Frame tick counter for LRU recycling */
	uint8_t        state;           /* 0 = free, 1 = ready */
	uint32_t       flags;           /* Aggregated node bitmask */
	int16_t        min_val;         /* Minimum scalar value in node */
	int16_t        max_val;         /* Maximum scalar value in node */
	uint32_t       rev;             /* Content revision counter */
	uint64_t       morton_key;      /* Multi-scale Morton key from morton_quad_key */
	uint8_t        active_view_mask; /* Bitmask of views active this frame (1<<view_id) */
	uint32_t       last_touch_frame; /* Global frame counter when last evaluated */
	uint16_t       grid_dim;        /* Runtime grid_dim for this node */
	uint16_t       vert_dim;        /* Runtime vert_dim (grid_dim + 1) */
	uint32_t       vert_count;      /* Runtime vert_count (vert_dim * vert_dim) */

	/* Pointers into the single contiguous arena allocated by ar_system_create_ex */
	int16_t       *values;          /* [vert_count] scalar values for zero-copy GPU upload */
	uint32_t      *rgba;            /* [vert_count] packed RGBA values */
	ar_sample_t   *samples;         /* [vert_count] full vertex samples */
} ar_node_t;

/* 6-plane frustum (normalized planes: A*x + B*y + C*z + D >= 0 is inside) */
typedef struct {
	float planes[6][4]; /* Normalized [A, B, C, D] for Left, Right, Bottom, Top, Near, Far */
} ar_frustum_t;

/* Options for matrix frustum culling, vertical envelope, and LOD splitting */
typedef struct {
	float        z_min, z_max;      /* Vertical bounds envelope */
	float        sag_radius;        /* Optional spherical sag radius R (0 = flat or use bounds_fn) */
	float        focal_x, focal_y;  /* Focal center for distance priority and sag computation */
	float        margin;            /* World-space padding margin added to node bounds */
	float        target_quad_px;    /* Screen-space quad LOD subdivision threshold in pixels */
	float        viewport_w;        /* Viewport width in pixels */
	float        viewport_h;        /* Viewport height in pixels */
	uint8_t      min_level;         /* Minimum LOD level clamp */
	uint8_t      max_level;         /* Maximum LOD level clamp */
	uint8_t      view_id;           /* View slot identifier (0..7) for multi-view active mask */
	uint16_t     max_populates_per_query; /* Optional per-query new-chunk populate cap (0 = unlimited) */
	float        split_hysteresis;  /* Optional split/merge hysteresis band (e.g. 0.12f) to prevent movement flicker */
	ar_bounds_fn bounds_fn;         /* Optional generic vertical envelope callback */
	void        *bounds_ctx;        /* Opaque context passed to bounds_fn */
} ar_cull_opts_t;

/* Camera viewpoint and frustum descriptor */
typedef struct {
	float   cam_x, cam_y, cam_z;     /* Camera focal position in world units */
	float   pitch;                   /* Pitch in radians (0 = top-down, ~0.96 = 55 deg diorama) */
	float   yaw;                     /* Yaw in radians */
	float   fov;                     /* Field of view in radians (0 = orthographic) */
	float   zoom;                    /* Global zoom multiplier (>0.0f) */
	float   viewport_w, viewport_h;  /* Viewport dimensions in pixels */
	float   eff_cam_dist;            /* Effective eye-to-target distance */
	uint8_t min_level;               /* Minimum LOD level clamp */
	uint8_t max_level;               /* Maximum LOD level clamp */
} ar_view_t;

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

bool ar_system_probe_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                          int32_t wx, int32_t wy,
                          ar_sample_t *out);

const ar_node_t *ar_system_find(ar_system_t *sys, uint16_t world,
                                uint8_t level, int32_t ox, int32_t oy);

const ar_node_t *ar_system_find_by_key(ar_system_t *sys, uint64_t morton_key);

const ar_node_t *ar_system_ensure(ar_system_t *sys, uint16_t world,
                                  uint8_t level, int32_t ox, int32_t oy);

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

static inline float ar_elevation_to_float(int16_t q12_4)
{
	return (float)q12_4 / 16.0f;
}

static inline int16_t ar_elevation_from_float(float val)
{
	return (int16_t)(val * 16.0f);
}

/**
 * Compute the world-space sampling step (texel span) for a quadtree LOD level.
 * Used to pass into chaff_band_range_for_step_q16 for Nyquist frequency culling.
 */
static inline float ar_level_texel_span(int32_t base_span, uint16_t grid_dim, uint8_t level)
{
	if (grid_dim == 0 || level >= 31)
		return 0.0f;
	float chunk_span = (float)(base_span >> level);
	return chunk_span / (float)grid_dim;
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
