#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ttypt/areole.h>
#include <ttypt/iso.h>

#define AR_WASM_DEFAULT_GRID_DIM     32u
#define AR_WASM_DEFAULT_POOL_CAP     128u
#define AR_WASM_DEFAULT_MAX_ACTIVE   64u
#define AR_WASM_DEFAULT_MAX_LEVELS   8u
#define AR_WASM_DEFAULT_RGBA         0xFF808080u

/* Generic analytical continuous scalar field sampler for standalone WASM queries */
static void default_wasm_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                 int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	(void)world;
	(void)level;

	float fx = (float)wx;
	float fy = (float)wy;

	float coarse = 120.0f * sinf(fx * 0.0030679f) * cosf(fy * 0.0030679f);
	float medium =  25.0f * sinf(fx * 0.049087f)  * cosf(fy * 0.049087f);
	float fine   =   6.0f * sinf(fx * 0.19635f)   * cosf(fy * 0.19635f);
	float field  = coarse + medium + fine;

	memset(out, 0, sizeof(*out));
	out->val = ar_value_from_float(field);
	out->rgba = AR_WASM_DEFAULT_RGBA;
	out->flags = (field < 0.0f) ? AR_NODE_FLAG_USER0 : 0u;
}

ISO_EXPORT(ar_system_t *, ar_wasm_create, int32_t base_w, int32_t base_h,
           uint32_t grid_dim, uint32_t pool_capacity, uint32_t max_active, uint32_t max_levels) {
	if (!grid_dim) grid_dim = AR_WASM_DEFAULT_GRID_DIM;
	if (!pool_capacity) pool_capacity = AR_WASM_DEFAULT_POOL_CAP;
	if (!max_active) max_active = AR_WASM_DEFAULT_MAX_ACTIVE;
	if (!max_levels) max_levels = AR_WASM_DEFAULT_MAX_LEVELS;
	return ar_system_create(base_w, base_h,
	                        (uint16_t)grid_dim,
	                        (uint16_t)pool_capacity,
	                        (uint16_t)max_active,
	                        (uint8_t)max_levels,
	                        default_wasm_sampler, NULL, NULL);
}

ISO_EXPORT(void, ar_wasm_destroy, ar_system_t *sys) {
	if (sys) {
		ar_system_destroy(sys);
	}
}

ISO_EXPORT(void, ar_wasm_clear, ar_system_t *sys) {
	if (sys) ar_system_clear(sys);
}

ISO_EXPORT(uint32_t, ar_wasm_update, ar_system_t *sys, uint32_t world,
                          float cam_x, float cam_y, float cam_z,
                          float pitch, float yaw, float fov,
                          float zoom, float vp_w, float vp_h, float eff_dist) {
	if (!sys) return 0;
	ar_view_t view = {
		.cam_x = cam_x, .cam_y = cam_y, .cam_z = cam_z,
		.pitch = pitch, .yaw = yaw, .fov = fov,
		.zoom = zoom,
		.viewport_w = vp_w, .viewport_h = vp_h,
		.eff_cam_dist = eff_dist
	};
	ar_stats_t st;
	ar_system_get_stats(sys, &st);
	size_t cap = st.pool_capacity > 0 ? st.pool_capacity : AR_WASM_DEFAULT_MAX_ACTIVE;
	const ar_node_t **scratch = (const ar_node_t **)malloc(cap * sizeof(const ar_node_t *));
	if (!scratch) return 0;
	size_t count = ar_system_update(sys, (uint16_t)world, &view, scratch, cap);
	free(scratch);
	return (uint32_t)count;
}

ISO_EXPORT(uint32_t, ar_wasm_get_system_node_count, const ar_system_t *sys) {
	return (uint32_t)ar_system_get_active_count(sys);
}

ISO_EXPORT(const ar_node_t *, ar_wasm_get_system_node_ptr, const ar_system_t *sys, uint32_t index) {
	return ar_system_get_active_node(sys, (size_t)index);
}

ISO_EXPORT(uint32_t, ar_wasm_get_node_count, const ar_system_t *sys) {
	return (uint32_t)ar_system_get_active_count(sys);
}

ISO_EXPORT(const ar_node_t *, ar_wasm_get_node_ptr, const ar_system_t *sys, uint32_t index) {
	return ar_system_get_active_node(sys, (size_t)index);
}

ISO_EXPORT(int32_t, ar_wasm_probe, ar_system_t *sys, uint32_t world, int32_t wx, int32_t wy, ar_sample_t *out) {
	if (!sys || !out) return 0;
	return ar_system_probe(sys, (uint16_t)world, wx, wy, out) ? 1 : 0;
}
