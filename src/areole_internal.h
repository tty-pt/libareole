#ifndef AREOLE_INTERNAL_H
#define AREOLE_INTERNAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if !defined(__wasm__) && (defined(__unix__) || defined(__APPLE__))
#include <time.h>
#endif

#include "ttypt/morton.h"
#include "ttypt/islet_quad.h"
#include "ttypt/areole.h"

typedef struct {
	int32_t ox, oy;
	uint8_t level;
	float   proj_w;
	float   dist_sq;
	float   split_ratio;
	float   b_min[3];
	float   b_max[3];
	bool    active;
	bool    split_blocked;
	bool    cached;
} ar_candidate_t;

typedef struct {
	uint64_t       morton_key;
	uint32_t       last_used_tick;
	uint32_t       last_touch_frame;
	uint32_t       active_view_mask;
	uint32_t       prev_active_view_mask;
	uint8_t        state;
	int32_t        sample_slot;
	int32_t        lru_prev;
	int32_t        lru_next;
#ifdef TTYPT_WQ_H
	stl_wq_task_t  async_task;
	struct ar_system *async_sys;
	uint32_t       async_slot;
	uint64_t       async_key;
	uint16_t       async_world;
	uint8_t        async_level;
	int32_t        async_ox;
	int32_t        async_oy;
	int32_t        async_w;
	int32_t        async_h;
	uint32_t       async_flags;
	int16_t        async_min_val;
	int16_t        async_max_val;
	int16_t       *staging_values;
	uint32_t      *staging_rgba;
	ar_sample_t   *staging_samples;
#endif
} ar_node_meta_t;

struct ar_system {
	int32_t           base_w;
	int32_t           base_h;
	uint16_t          grid_dim;
	uint16_t          vert_dim;
	uint32_t          vert_count;
	uint16_t          pool_capacity;
	uint16_t          samples_pool_capacity;
	uint8_t           samples_max_level;
	bool              samples_bounded;
	uint16_t          max_active;
	uint8_t           max_levels;
	uint32_t          default_rgba;
	float             z_min;
	float             z_max;
	float             target_quad_px;
	float             margin;
	float             value_scale;
	uint16_t          max_root_cells;
	uint8_t           probe_ancestor_cap;
	int32_t           period_x;
	int32_t           period_y;
	ar_sampler_fn     sampler;
	ar_row_sampler_fn row_sampler;
	ar_sparsity_fn    sparsity;
	void             *user_ctx;
	uint32_t          max_populates_per_frame;
	uint32_t          frame_populates_count;
	uint32_t          map_capacity;
	ar_interp_mode_t  interp_mode;
	ar_interp_fn      interp_fn;
	float             default_pitch;
	float             default_cam_dist;
	float             default_topdown_dist;
	float             default_viewport_w;
	float             default_viewport_h;
	float             sag_radius;
	float             sag_tangent_slope;
	float             split_hysteresis;
	ar_bounds_fn      bounds_fn;
	ar_bounds_3d_fn   bounds_3d_fn;
	ar_cull_3d_fn     cull_3d_fn;
	void             *bounds_ctx;
	uint8_t           attr_count;
	ar_interp_mode_t  attr_interp_mode;
	bool              interpolate_attrs;
	uint32_t          current_tick;
	uint32_t          begin_frame_cleared_mask;
	struct stl_wq    *wq;

	/* O(1) free-list + intrusive doubly-linked LRU */
	int32_t           free_head;
	int32_t           lru_head;
	int32_t           lru_tail;
	size_t            pool_in_use;

	/* Pointers into single contiguous allocation */
	ar_node_t        *pool;
	ar_node_meta_t   *pool_meta;
	const ar_node_t **active_chunks;
	ar_candidate_t   *cands_scratch;
	ar_candidate_t   *active_cands_scratch;
	uint32_t         *heap_scratch;
	uint32_t         *free_scratch;
	uint8_t          *samples_arena;
	size_t            per_node_samples_stride;
	int32_t          *sample_owner;
	ar_sample_t      *row_scratch;
	size_t            active_count;
	uint8_t          *staging_arena;
	size_t            staging_bytes;

	/* Backing libislet in-memory quadtree spatial database */
	islet_quad_map_t  map;

	/* Telemetry & diagnostics */
	ar_stats_t        stats;
	ar_metrics_t      metrics;
	bool              timing_enabled;
};

static inline size_t ar_align8(size_t n)
{
	return (n + 7u) & ~(size_t)7u;
}

static inline uint64_t ar_now_us(bool enabled)
{
	if (!enabled)
		return 0u;
#if !defined(__wasm__) && (defined(__unix__) || defined(__APPLE__))
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (uint64_t)ts.tv_sec * 1000000ull + (uint64_t)ts.tv_nsec / 1000ull;
#else
	return 0u;
#endif
}

static inline bool ar_map_get(ar_system_t *sys, uint64_t key, uint32_t *out_slot)
{
	uint32_t steps = 0;
	bool found = islet_quad_map_get_probed(&sys->map, key, out_slot, &steps);
	sys->metrics.frame_map_lookups++;
	sys->metrics.frame_map_probe_steps += steps;
	return found;
}

/* Internal cross-unit declarations */
static void ar_node_compute_occ_cell_min(ar_node_t *node);
static void ar_node_release_samples(ar_system_t *sys, int slot);
static void ar_node_bind_samples(ar_system_t *sys, int slot, uint8_t level);
static void ar_lru_touch(ar_system_t *sys, int32_t slot);
static void ar_lru_insert_head(ar_system_t *sys, int32_t slot);
static int  ar_pool_alloc_or_evict(ar_system_t *sys);
static ar_node_t *ar_get_or_populate(ar_system_t *sys, uint16_t world,
                                     uint8_t level, int32_t ox, int32_t oy);
#ifdef TTYPT_WQ_H
static void ar_async_populate_worker(void *ctx);
#endif

#endif /* AREOLE_INTERNAL_H */
