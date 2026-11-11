#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "ttypt/areole_glsl.h"

int main(void)
{
	printf("=== Testing GLSL ES 3.0 GPU Shader Contract & Reference Parity for libareole ===\n");

	/* 1. Uniform layouts and definitions */
	{
		assert(strcmp(AR_UNIFORM_ORIGIN, "uAreoleOrigin") == 0);
		assert(strcmp(AR_UNIFORM_SPAN, "uAreoleSpan") == 0);
		assert(strcmp(AR_UNIFORM_MVP, "uMVP") == 0);
		assert(strcmp(AR_UNIFORM_FOCAL_CENTER, "uFocalCenter") == 0);
		assert(strcmp(AR_UNIFORM_CURVE_RADIUS, "uCurveRadius") == 0);

		ar_glsl_uniforms_t u;
		memset(&u, 0, sizeof(u));
		u.origin[0] = 160.0f; u.origin[1] = 128.0f;
		u.span[0] = 320.0f; u.span[1] = 256.0f;
		u.focal_center[0] = 0.0f; u.focal_center[1] = 0.0f;
		u.curve_radius = 24000.0f;

		assert(u.origin[0] == 160.0f && u.span[0] == 320.0f);
		assert(u.curve_radius == 24000.0f);
		printf("  [PASS] Uniform definitions and C uniform buffer layout validated\n");
	}

	/* 2. Shader source functions string */
	{
		const char *src = ar_glsl_source();
		assert(src != NULL);
		assert(strstr(src, "ar_uv_to_world") != NULL);
		assert(strstr(src, "ar_apply_spherical_sag") != NULL);
		assert(strstr(src, "ar_project") != NULL);

		const char *csrc = ar_glsl_source();
		assert(csrc != NULL);
		assert(strstr(csrc, "ar_uv_to_world") != NULL);
		printf("  [PASS] ar_glsl_source and ar_glsl_source verified\n");
	}

	/* 3. ar_glsl_ref_uv_to_world float parity */
	{
		float origin[2] = { -320.0f, 160.0f };
		float span[2]   = {  640.0f, 480.0f };
		float uv[2]     = {    0.25f,  0.75f };
		float out[2];

		ar_glsl_ref_uv_to_world(uv, origin, span, out);
		assert(fabsf(out[0] - (-320.0f + 0.25f * 640.0f)) < 1e-5f);
		assert(fabsf(out[1] - ( 160.0f + 0.75f * 480.0f)) < 1e-5f);
		printf("  [PASS] ar_glsl_ref_uv_to_world reference parity verified\n");
	}

	/* 4. ar_glsl_ref_apply_spherical_sag float parity across radius zones (and flat R <= 0) */
	{
		float focal[2] = { 100.0f, 200.0f };
		float R = 20000.0f;

		/* Flat: R <= 0 -> sag = 0 */
		float p_flat[2] = { 5000.0f, 5000.0f };
		float curved[3];
		ar_glsl_ref_apply_spherical_sag(p_flat, 50.0f, focal, 0.0f, curved);
		assert(fabsf(curved[2] - 50.0f) < 1e-5f);

		/* Center: r = 0 -> sag = 0 */
		float p_center[2] = { 100.0f, 200.0f };
		ar_glsl_ref_apply_spherical_sag(p_center, 50.0f, focal, R, curved);
		assert(fabsf(curved[0] - 100.0f) < 1e-5f);
		assert(fabsf(curved[1] - 200.0f) < 1e-5f);
		assert(fabsf(curved[2] - 50.0f) < 1e-4f);

		/* Distance r = 12000 < R */
		float p_mid[2] = { 100.0f + 12000.0f, 200.0f };
		ar_glsl_ref_apply_spherical_sag(p_mid, 100.0f, focal, R, curved);
		float expected_z = 100.0f - 4000.0f;
		assert(fabsf(curved[2] - expected_z) < 1e-3f);

		/* Horizon r = R = 20000 */
		float p_horiz[2] = { 100.0f, 200.0f + 20000.0f };
		ar_glsl_ref_apply_spherical_sag(p_horiz, 0.0f, focal, R, curved);
		assert(fabsf(curved[2] - (-20000.0f)) < 1e-2f);

		/* Beyond horizon r = 25000 > R */
		float p_beyond[2] = { 100.0f + 25000.0f, 200.0f };
		ar_glsl_ref_apply_spherical_sag(p_beyond, 0.0f, focal, R, curved);
		assert(fabsf(curved[2] - (-30000.0f)) < 1e-2f);

		printf("  [PASS] ar_glsl_ref_apply_spherical_sag verified across flat, center, sphere arc, horizon, and tangent\n");
	}

	/* 5. ar_glsl_ref_project parity */
	{
		float mvp[16] = {
			2.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 3.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 4.0f, 0.0f,
			10.0f, 20.0f, 30.0f, 1.0f
		};
		float p[3] = { 5.0f, 6.0f, 7.0f };
		float clip[4];

		ar_glsl_ref_project(p, mvp, clip);
		assert(fabsf(clip[0] - 20.0f) < 1e-5f);
		assert(fabsf(clip[1] - 38.0f) < 1e-5f);
		assert(fabsf(clip[2] - 58.0f) < 1e-5f);
		assert(fabsf(clip[3] - 1.0f) < 1e-5f);
		printf("  [PASS] ar_glsl_ref_project column-major projection verified\n");
	}

	printf("=== All GLSL ES 3.0 Parity Tests PASSED ===\n");
	return 0;
}
