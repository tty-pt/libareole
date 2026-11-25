/* ========================================================================= */
/* Matrix Mathematics, Frustum Extraction & Projection Helpers               */
/* ========================================================================= */

void ar_mat4_identity(float out[16])
{
	memset(out, 0, 16 * sizeof(float));
	out[0] = 1.0f;
	out[5] = 1.0f;
	out[10] = 1.0f;
	out[15] = 1.0f;
}

void ar_mat4_multiply(float out[16], const float a[16], const float b[16])
{
	float tmp[16];
	for (int c = 0; c < 4; c++) {
		for (int r = 0; r < 4; r++) {
			tmp[c * 4 + r] = a[0 * 4 + r] * b[c * 4 + 0] +
			                 a[1 * 4 + r] * b[c * 4 + 1] +
			                 a[2 * 4 + r] * b[c * 4 + 2] +
			                 a[3 * 4 + r] * b[c * 4 + 3];
		}
	}
	memcpy(out, tmp, 16 * sizeof(float));
}

void ar_mat4_ortho(float out[16], float left, float right, float bottom, float top, float near_z, float far_z)
{
	memset(out, 0, 16 * sizeof(float));
	float rl = right - left;
	float tb = top - bottom;
	float fn = far_z - near_z;
	if (fabsf(rl) < 1e-7f || fabsf(tb) < 1e-7f || fabsf(fn) < 1e-7f) {
		ar_mat4_identity(out);
		return;
	}
	out[0]  = 2.0f / rl;
	out[5]  = 2.0f / tb;
	out[10] = -2.0f / fn;
	out[12] = -(right + left) / rl;
	out[13] = -(top + bottom) / tb;
	out[14] = -(far_z + near_z) / fn;
	out[15] = 1.0f;
}

void ar_mat4_perspective(float out[16], float fov_rad, float aspect, float near_z, float far_z)
{
	memset(out, 0, 16 * sizeof(float));
	float half_fov = fov_rad * 0.5f;
	float sin_fov = sinf(half_fov);
	float cos_fov = cosf(half_fov);
	if (fabsf(sin_fov) < 1e-7f || fabsf(aspect) < 1e-7f || fabsf(far_z - near_z) < 1e-7f) {
		ar_mat4_identity(out);
		return;
	}
	float f = cos_fov / sin_fov;
	float nf = 1.0f / (near_z - far_z);

	out[0]  = f / aspect;
	out[5]  = f;
	out[10] = (far_z + near_z) * nf;
	out[11] = -1.0f;
	out[14] = 2.0f * far_z * near_z * nf;
}

/* Standard right-handed OpenGL/WebGL LookAt */
void ar_mat4_lookat(float out[16], const float eye[3], const float target[3], const float up[3])
{
	float z0 = eye[0] - target[0];
	float z1 = eye[1] - target[1];
	float z2 = eye[2] - target[2];
	float len_z = sqrtf(z0 * z0 + z1 * z1 + z2 * z2);
	if (len_z < 1e-7f) {
		z0 = 0.0f; z1 = 0.0f; z2 = 1.0f;
	} else {
		float inv_z = 1.0f / len_z;
		z0 *= inv_z; z1 *= inv_z; z2 *= inv_z;
	}

	float upx = up[0], upy = up[1], upz = up[2];
	float x0 = upy * z2 - upz * z1;
	float x1 = upz * z0 - upx * z2;
	float x2 = upx * z1 - upy * z0;
	float len_x = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);
	if (len_x < 1e-7f) {
		x0 = 1.0f; x1 = 0.0f; x2 = 0.0f;
	} else {
		float inv_x = 1.0f / len_x;
		x0 *= inv_x; x1 *= inv_x; x2 *= inv_x;
	}

	float y0 = z1 * x2 - z2 * x1;
	float y1 = z2 * x0 - z0 * x2;
	float y2 = z0 * x1 - z1 * x0;

	out[0] = x0; out[1] = y0; out[2] = z0; out[3] = 0.0f;
	out[4] = x1; out[5] = y1; out[6] = z1; out[7] = 0.0f;
	out[8] = x2; out[9] = y2; out[10] = z2; out[11] = 0.0f;
	out[12] = -(x0 * eye[0] + x1 * eye[1] + x2 * eye[2]);
	out[13] = -(y0 * eye[0] + y1 * eye[1] + y2 * eye[2]);
	out[14] = -(z0 * eye[0] + z1 * eye[1] + z2 * eye[2]);
	out[15] = 1.0f;
}

/* Left-handed / East-flipped LookAt convention */
void ar_mat4_lookat_lh(float out[16], const float eye[3], const float target[3], const float up[3])
{
	float z0 = eye[0] - target[0];
	float z1 = eye[1] - target[1];
	float z2 = eye[2] - target[2];
	float len_z = sqrtf(z0 * z0 + z1 * z1 + z2 * z2);
	if (len_z < 1e-7f) {
		z0 = 0.0f; z1 = 0.0f; z2 = 1.0f;
	} else {
		float inv_z = 1.0f / len_z;
		z0 *= inv_z; z1 *= inv_z; z2 *= inv_z;
	}

	float upx = up[0], upy = up[1], upz = up[2];
	float x0 = -(upy * z2 - upz * z1);
	float x1 = -(upz * z0 - upx * z2);
	float x2 = -(upx * z1 - upy * z0);
	float len_x = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);
	if (len_x < 1e-7f) {
		x0 = 1.0f; x1 = 0.0f; x2 = 0.0f;
	} else {
		float inv_x = 1.0f / len_x;
		x0 *= inv_x; x1 *= inv_x; x2 *= inv_x;
	}

	float y0 = x1 * z2 - x2 * z1;
	float y1 = x2 * z0 - x0 * z2;
	float y2 = x0 * z1 - x1 * z0;

	out[0] = x0; out[1] = y0; out[2] = z0; out[3] = 0.0f;
	out[4] = x1; out[5] = y1; out[6] = z1; out[7] = 0.0f;
	out[8] = x2; out[9] = y2; out[10] = z2; out[11] = 0.0f;
	out[12] = -(x0 * eye[0] + x1 * eye[1] + x2 * eye[2]);
	out[13] = -(y0 * eye[0] + y1 * eye[1] + y2 * eye[2]);
	out[14] = -(z0 * eye[0] + z1 * eye[1] + z2 * eye[2]);
	out[15] = 1.0f;
}

/* 4x4 Matrix Inversion (double-precision cofactors for high near/far ratios) */
bool ar_mat4_invert(float out[16], const float mf[16])
{
	double m[16];
	for (int i = 0; i < 16; i++)
		m[i] = (double)mf[i];

	double inv[16];
	inv[0] = m[5]  * m[10] * m[15] -
	         m[5]  * m[11] * m[14] -
	         m[9]  * m[6]  * m[15] +
	         m[9]  * m[7]  * m[14] +
	         m[13] * m[6]  * m[11] -
	         m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] +
	          m[4]  * m[11] * m[14] +
	          m[8]  * m[6]  * m[15] -
	          m[8]  * m[7]  * m[14] -
	          m[12] * m[6]  * m[11] +
	          m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] -
	         m[4]  * m[11] * m[13] -
	         m[8]  * m[5] * m[15] +
	         m[8]  * m[7] * m[13] +
	         m[12] * m[5] * m[11] -
	         m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] +
	           m[4]  * m[10] * m[13] +
	           m[8]  * m[5] * m[14] -
	           m[8]  * m[6] * m[13] -
	           m[12] * m[5] * m[10] +
	           m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] +
	          m[1]  * m[11] * m[14] +
	          m[9]  * m[2] * m[15] -
	          m[9]  * m[3] * m[14] -
	          m[13] * m[2] * m[11] +
	          m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] -
	         m[0]  * m[11] * m[14] -
	         m[8]  * m[2] * m[15] +
	         m[8]  * m[3] * m[14] +
	         m[12] * m[2] * m[11] -
	         m[12] * m[3] * m[10];

	inv[9] = -m[0]  * m[9] * m[15] +
	          m[0]  * m[11] * m[13] +
	          m[8]  * m[1] * m[15] -
	          m[8]  * m[3] * m[13] -
	          m[12] * m[1] * m[11] +
	          m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] -
	          m[0]  * m[10] * m[13] -
	          m[8]  * m[1] * m[14] +
	          m[8]  * m[2] * m[13] +
	          m[12] * m[1] * m[10] -
	          m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] -
	         m[1]  * m[7] * m[14] -
	         m[5]  * m[2] * m[15] +
	         m[5]  * m[3] * m[14] +
	         m[13] * m[2] * m[7] -
	         m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] +
	          m[0]  * m[7] * m[14] +
	          m[4]  * m[2] * m[15] -
	          m[4]  * m[3] * m[14] -
	          m[12] * m[2] * m[7] +
	          m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] -
	          m[0]  * m[7] * m[13] -
	          m[4]  * m[1] * m[15] +
	          m[4]  * m[3] * m[13] +
	          m[12] * m[1] * m[7] -
	          m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] +
	           m[0]  * m[6] * m[13] +
	           m[4]  * m[1] * m[14] -
	           m[4]  * m[2] * m[13] -
	           m[12] * m[1] * m[6] +
	           m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
	          m[1] * m[7] * m[10] +
	          m[5] * m[2] * m[11] -
	          m[5] * m[3] * m[10] -
	          m[9] * m[2] * m[7] +
	          m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
	         m[0] * m[7] * m[10] -
	         m[4] * m[2] * m[11] +
	         m[4] * m[3] * m[10] +
	         m[8] * m[2] * m[7] -
	         m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
	           m[0] * m[7] * m[9] +
	           m[4] * m[1] * m[11] -
	           m[4] * m[3] * m[9] -
	           m[8] * m[1] * m[7] +
	           m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
	          m[0] * m[6] * m[9] -
	          m[4] * m[1] * m[10] +
	          m[4] * m[2] * m[9] +
	          m[8] * m[1] * m[6] -
	          m[8] * m[2] * m[5];

	double det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
	if (fabs(det) < 1e-35 || !isfinite(det)) {
		ar_mat4_identity(out);
		return false;
	}

	double inv_det = 1.0 / det;
	for (int i = 0; i < 16; i++) {
		out[i] = (float)(inv[i] * inv_det);
	}
	return true;
}

void ar_frustum_from_matrix(ar_frustum_t *out_frustum, const float m[16])
{
	if (!out_frustum || !m) return;

	/* Left: row3 + row0 */
	out_frustum->planes[0][0] = m[3]  + m[0];
	out_frustum->planes[0][1] = m[7]  + m[4];
	out_frustum->planes[0][2] = m[11] + m[8];
	out_frustum->planes[0][3] = m[15] + m[12];

	/* Right: row3 - row0 */
	out_frustum->planes[1][0] = m[3]  - m[0];
	out_frustum->planes[1][1] = m[7]  - m[4];
	out_frustum->planes[1][2] = m[11] - m[8];
	out_frustum->planes[1][3] = m[15] - m[12];

	/* Bottom: row3 + row1 */
	out_frustum->planes[2][0] = m[3]  + m[1];
	out_frustum->planes[2][1] = m[7]  + m[5];
	out_frustum->planes[2][2] = m[11] + m[9];
	out_frustum->planes[2][3] = m[15] + m[13];

	/* Top: row3 - row1 */
	out_frustum->planes[3][0] = m[3]  - m[1];
	out_frustum->planes[3][1] = m[7]  - m[5];
	out_frustum->planes[3][2] = m[11] - m[9];
	out_frustum->planes[3][3] = m[15] - m[13];

	/* Near: row3 + row2 (OpenGL NDC z in [-1, 1]) */
	out_frustum->planes[4][0] = m[3]  + m[2];
	out_frustum->planes[4][1] = m[7]  + m[6];
	out_frustum->planes[4][2] = m[11] + m[10];
	out_frustum->planes[4][3] = m[15] + m[14];

	/* Far: row3 - row2 */
	out_frustum->planes[5][0] = m[3]  - m[2];
	out_frustum->planes[5][1] = m[7]  - m[6];
	out_frustum->planes[5][2] = m[11] - m[10];
	out_frustum->planes[5][3] = m[15] - m[14];

	for (int i = 0; i < 6; i++) {
		float a = out_frustum->planes[i][0];
		float b = out_frustum->planes[i][1];
		float c = out_frustum->planes[i][2];
		float len = sqrtf(a * a + b * b + c * c);
		if (len > 1e-7f) {
			float inv_len = 1.0f / len;
			out_frustum->planes[i][0] *= inv_len;
			out_frustum->planes[i][1] *= inv_len;
			out_frustum->planes[i][2] *= inv_len;
			out_frustum->planes[i][3] *= inv_len;
		}
	}
}

bool ar_frustum_test_box(const ar_frustum_t *frustum,
                         float min_x, float min_y, float min_z,
                         float max_x, float max_y, float max_z)
{
	if (!frustum) return true;

	for (int i = 0; i < 6; i++) {
		float a = frustum->planes[i][0];
		float b = frustum->planes[i][1];
		float c = frustum->planes[i][2];
		float d = frustum->planes[i][3];

		float px = (a > 0.0f) ? max_x : min_x;
		float py = (b > 0.0f) ? max_y : min_y;
		float pz = (c > 0.0f) ? max_z : min_z;

		if (a * px + b * py + c * pz + d < 0.0f) {
			return false;
		}
	}
	return true;
}

static inline float ar_projected_screen_width_box_3d(const float mvp[16],
                                                     const ar_cull_opts_t *opts,
                                                     const float center_xyz[3],
                                                     const float eye_xyz[3], bool has_eye,
                                                     float span_w, float span_h,
                                                     float vp_w, float vp_h)
{
	bool is_perspective = (fabsf(mvp[3]) > 1e-5f || fabsf(mvp[7]) > 1e-5f || fabsf(mvp[11]) > 1e-5f);

	float len_row0 = sqrtf(mvp[0]*mvp[0] + mvp[4]*mvp[4] + mvp[8]*mvp[8]);
	float len_row1 = sqrtf(mvp[1]*mvp[1] + mvp[5]*mvp[5] + mvp[9]*mvp[9]);
	if (len_row0 < 1e-7f) len_row0 = 0.001f;
	if (len_row1 < 1e-7f) len_row1 = 0.001f;

	float w_eff = 1.0f;
	if (is_perspective) {
		float cw = mvp[3] * center_xyz[0] +
		           mvp[7] * center_xyz[1] +
		           mvp[11] * center_xyz[2] +
		           mvp[15];
		float clip_w_thresh = (opts && opts->behind_cam_clip_w > 0.0f)
		                      ? opts->behind_cam_clip_w : AR_DEFAULT_BEHIND_CAM_CLIP_W;
		if (has_eye) {
			float dx = center_xyz[0] - eye_xyz[0];
			float dy = center_xyz[1] - eye_xyz[1];
			float dz = center_xyz[2] - eye_xyz[2];
			float len_row3 = sqrtf(mvp[3]*mvp[3] + mvp[7]*mvp[7] + mvp[11]*mvp[11]);
			float d_center = sqrtf(dx * dx + dy * dy + dz * dz) * len_row3;
			w_eff = (d_center > clip_w_thresh) ? d_center : clip_w_thresh;
		} else {
			w_eff = (cw > clip_w_thresh) ? cw : clip_w_thresh;
		}
	}

	float px_w = (span_w * len_row0 * vp_w * 0.5f) / w_eff;
	float px_h = (span_h * len_row1 * vp_h * 0.5f) / w_eff;

	return (px_w > px_h) ? px_w : px_h;
}

static inline bool ar_unproject_ndc(const float inv_mvp[16], float ndc_x, float ndc_y, float ndc_z, float out_world[3])
{
	float v[4] = { ndc_x, ndc_y, ndc_z, 1.0f };
	float p[4];
	for (int i = 0; i < 4; i++) {
		p[i] = inv_mvp[0 * 4 + i] * v[0] +
		       inv_mvp[1 * 4 + i] * v[1] +
		       inv_mvp[2 * 4 + i] * v[2] +
		       inv_mvp[3 * 4 + i] * v[3];
	}
	if (fabsf(p[3]) < 1e-9f) return false;
	float inv_w = 1.0f / p[3];
	out_world[0] = p[0] * inv_w;
	out_world[1] = p[1] * inv_w;
	out_world[2] = p[2] * inv_w;
	return true;
}

static inline int32_t ar_floor_i32(float v)
{
	if (v < -1e9f) return -1000000000;
	if (v >  1e9f) return  1000000000;
	return (int32_t)floorf(v);
}

static inline bool ar_frustum_test_box_sat(const ar_frustum_t *f,
                                           bool has_faabb,
                                           const float f_min[3], const float f_max[3],
                                           const float b_min[3], const float b_max[3])
{
	if (has_faabb) {
		if (b_max[0] < f_min[0] || b_min[0] > f_max[0] ||
		    b_max[1] < f_min[1] || b_min[1] > f_max[1] ||
		    b_max[2] < f_min[2] || b_min[2] > f_max[2])
			return false;
	}
	return ar_frustum_test_box(f, b_min[0], b_min[1], b_min[2],
	                           b_max[0], b_max[1], b_max[2]);
}
