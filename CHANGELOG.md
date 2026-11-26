# Changelog

## [Unreleased] (target: v1.2.0)
- **4-Wide Batched Row Populate Callback (`ar_eval_row4_fn`, `include/ttypt/areole.h`, `src/ar_populate.c`)**:
  - Added optional `ar_eval_row4_fn` (`cfg.eval_row4`) for 4-wide row-coherent batch sampling on mixed-reuse rows during chunk population, with parity tests in `tests/test_areole.c`.
- **6-Unit Modular Architecture, O(1) Pool/LRU, Telemetry & Breadth-First Budget Coarsening (`src/ar_*.c`, `src/areole_internal.h`, `include/ttypt/areole.h`, `tests/bench_areole.c`)**:
  - Decoupled monolithic `src/areole.c` into 6 focused compilation units: `ar_math.c`, `ar_occ.c`, `ar_pool.c`, `ar_populate.c`, `ar_probe.c`, and `ar_solver.c` (with shared `src/areole_internal.h` and `tests/bench_areole.c` benchmark harness).
  - Replaced linear pool scans with O(1) pool/LRU management, binary max-heap LOD refinement, and breadth-first budget coarsening.
  - Added detailed per-frame and structural telemetry (`ar_metrics_t`, `ar_system_get_metrics()`, `ar_system_reset_frame_metrics()`, `ar_system_set_timing_enabled()`).
- **Bounded `ar_sample_t` Sub-Pool & Arena Telemetry (`include/ttypt/areole.h`, `src/ar_pool.c`)**:
  - Added `samples_pool_capacity` and `samples_max_level` to `ar_config_t` to cap nodes retaining full `ar_sample_t[]` arrays (cutting peak arena footprint by ~20.5 MB on high-LOD configs) and added `arena_bytes` reporting to `ar_stats_t`.
- **Exact SAT Frustum Culling, Cached Z Tightening, 2D Software Hi-Z Occluder & Height-Only Probes (`include/ttypt/areole.h`, `src/ar_occ.c`, `src/ar_probe.c`, `src/ar_solver.c`)**:
  - Added exact Separating Axis Theorem (SAT) frustum culling, optional cached-node vertical bounds tightening (`use_cached_bounds`, `cached_bounds_max_level`, `cached_bounds_pad`), and `focal_z` camera ground-plane intersection in `ar_cull_opts_t`.
  - Added conservative 96×54 (8×8 tiled) 2D software Hi-Z linear-W occluder (`ar_occluder_t`, `ar_point_3d_fn`, `ar_occluder_reset()`, `ar_occluder_rasterize_tri()`, `ar_occluder_rasterize_node()`, `ar_occluder_test_box()`) with precomputed 4×4 `occ_cell_min` per node (`AR_OCC_CELL_K`).
  - Added height-only probe fast paths (`ar_system_probe_val()`, `ar_system_probe_view_val()`, `ar_probe_val()`, `ar_probe_view_val()`).
- **Scanline Row-Batch Sampler & Per-Frame Populate Budgeting (`ar_row_sampler_fn`, `include/ttypt/areole.h`, `src/ar_populate.c`)**:
  - Added optional `ar_row_sampler_fn` (`cfg.row_sampler`) scanline callback passing exact `span_w` world coordinates for vectorized row sampling and watertight cross-LOD boundary vertex welding.
  - Added `max_populates_per_frame` in `ar_config_t` and `max_populates_per_query` in `ar_cull_opts_t` to cap new-chunk synthesis per frame while falling back to cached ancestor LOD nodes.
- **Pluggable 3D Manifold Normal-Cone & Horizon Culling Callback (`ar_cull_3d_fn`, `include/ttypt/areole.h`, `src/ar_solver.c`)**:
  - Added `ar_cull_3d_fn` (`cull_3d_fn` in `ar_config_t` and `ar_cull_opts_t`) to cull back-facing or beyond-the-horizon curved manifold patches before quadtree subdivision and vertex population.
- **Asynchronous Workqueue (`stele_wq`) Chunk Population & Deferred Upsampling Fixes (`include/ttypt/areole.h`, `src/ar_populate.c`, `src/ar_solver.c`)**:
  - Integrated `<ttypt/wq.h>` (`struct stl_wq *wq` in `ar_config_t`, `AR_NODE_FLAG_DEFERRED`) so unpopulated leaf nodes submit zero-allocation background tasks (`stl_wq_submit`) and render parent ancestor fallbacks until worker population completes.
  - Fixed LOD population budget accounting and required populated parent nodes before deferred upsampling.

## [1.1.0] - 2026-11-16
- **Pluggable 3D Manifold Bounding Box Callback (`ar_bounds_3d_fn`, `include/ttypt/areole.h`, `src/areole.c`)**:
  - Added `ar_bounds_3d_fn` to `ar_config_t` and `ar_cull_opts_t` to map 2D manifold patches `[ox, ox + w] x [oy, oy + h]` and normal displacement `[z_min, z_max]` into true 3D Cartesian AABBs, centers, and chord spans for non-Euclidean and curved planetary frustum culling and LOD splitting.
- **Node Pool Encapsulation & Parameterization (`include/ttypt/areole.h`, `include/ttypt/areole_glsl.h`, `src/areole.c`)**:
  - Encapsulated internal pool metadata (`morton_key`, `active_views`, `in_use`) behind public accessors (`ar_system_node_morton_key`, `ar_system_node_active_views`, `ar_system_node_is_active_in_view`) and reserved `AR_NODE_FLAG_SPARSE` (`1u << 31`).
  - Parameterized probe interpolation (`ar_interp_mode_t`, `ar_interp_fn`, `attr_interp_mode`, `attr_count`), value scaling (`value_scale`, `ar_system_value_to_float`), camera/ray unprojection options, and GLSL spherical sag slope (`ar_apply_spherical_sag_ex`, `ar_glsl_ref_eval_vertex`).
  - Handled infinite far-plane unprojection in `ar_unproject_ndc` and tilted-frustum northern row intersection in `ar_system_query_matrix`.
- **Universal 7-OS & NPM CI Pipeline (`.github/workflows/ci.yml`, `Makefile`, `npm/native/areole_napi.c`)**:
  - Migrated `Makefile` to `include ../mk/include.mk`, added `npm/native/areole_napi.c`, and added `.github/workflows/ci.yml` publishing `libareole` across `deb,apk,rpm,pacman,brew,openbsd,winget` and `@tty-pt/areole` to npm with `linux-x64`, `darwin-arm64`, `darwin-x64`, `win32-x64` N-API prebuilds and `areole.wasm`.
- Linked `areole.wasm` / `chunk.wasm` with shared memory and the `@tty-pt/stele` single-allocator substrate (`>=1.7.0` / `^1.7.0`).

## [1.0.0] - 2026-04-19
- Initial public release of `libareole` (`<ttypt/areole.h>`, `<ttypt/areole_glsl.h>`, `@tty-pt/areole`).
- Runtime-configurable `ar_config_t` (`grid_dim`, `pool_capacity`, `max_active`, `max_levels`, `default_rgba`) with single contiguous arena allocation (`ar_system_create_ex`).
- Standard right-handed OpenGL/WebGL 4x4 matrix math (`ar_mat4_lookat`) and Gribb-Hartmann 6-plane frustum extraction (`ar_frustum_from_matrix`).
- Pluggable vertical envelope callback (`ar_bounds_fn`).
