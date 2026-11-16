# Changelog

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
