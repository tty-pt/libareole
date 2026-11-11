# Changelog

## [1.0.0] - 2026-04-19
- Initial public release of `libareole` (`<ttypt/areole.h>`, `<ttypt/areole_glsl.h>`, `@tty-pt/areole`).
- Runtime-configurable `ar_config_t` (`grid_dim`, `pool_capacity`, `max_active`, `max_levels`, `default_rgba`) with single contiguous arena allocation (`ar_system_create_ex`).
- Standard right-handed OpenGL/WebGL 4x4 matrix math (`ar_mat4_lookat`) and Gribb-Hartmann 6-plane frustum extraction (`ar_frustum_from_matrix`).
- Pluggable vertical envelope callback (`ar_bounds_fn`).
