#include <node_api.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ttypt/areole.h>

static void napi_default_sampler(void *user_ctx, uint16_t world, uint8_t level,
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
	out->rgba = 0xFF808080u;
	out->flags = (field < 0.0f) ? AR_NODE_FLAG_USER0 : 0u;
}

static napi_value AreoleMat4Identity(napi_env env, napi_callback_info info)
{
	(void)info;
	float m[16];
	ar_mat4_identity(m);
	void *data = NULL;
	napi_value ab, ta;
	napi_create_arraybuffer(env, sizeof(m), &data, &ab);
	if (data)
		memcpy(data, m, sizeof(m));
	napi_create_typedarray(env, napi_float32_array, 16, ab, 0, &ta);
	return ta;
}

static napi_value AreoleProbeDefault(napi_env env, napi_callback_info info)
{
	size_t argc = 3;
	napi_value args[3];
	napi_get_cb_info(env, info, &argc, args, NULL, NULL);

	int32_t wx = 0, wy = 0;
	uint32_t world = 1;
	if (argc >= 1) napi_get_value_int32(env, args[0], &wx);
	if (argc >= 2) napi_get_value_int32(env, args[1], &wy);
	if (argc >= 3) napi_get_value_uint32(env, args[2], &world);

	ar_system_t *sys = ar_system_create(160, 128, 32, 16, 8, 8, napi_default_sampler, NULL, NULL);
	if (!sys) return NULL;
	ar_sample_t s;
	memset(&s, 0, sizeof(s));
	bool ok = ar_system_probe(sys, (uint16_t)world, wx, wy, &s);
	ar_system_destroy(sys);
	if (!ok) return NULL;

	napi_value obj, v_val, v_rgba, v_flags;
	napi_create_object(env, &obj);
	napi_create_double(env, (double)s.val / 16.0, &v_val);
	napi_create_uint32(env, s.rgba, &v_rgba);
	napi_create_uint32(env, s.flags, &v_flags);
	napi_set_named_property(env, obj, "value", v_val);
	napi_set_named_property(env, obj, "rgba", v_rgba);
	napi_set_named_property(env, obj, "flags", v_flags);
	return obj;
}

static napi_value Init(napi_env env, napi_value exports)
{
	napi_property_descriptor desc[] = {
		{ "mat4Identity", NULL, AreoleMat4Identity, NULL, NULL, NULL, napi_default, NULL },
		{ "probeDefault", NULL, AreoleProbeDefault, NULL, NULL, NULL, napi_default, NULL }
	};
	napi_define_properties(env, exports, sizeof(desc) / sizeof(desc[0]), desc);
	return exports;
}

NAPI_MODULE(NODE_GYP_MODULE_NAME, Init)
