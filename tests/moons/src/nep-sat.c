#include <math.h>

static inline long double get_dfdx_cell(int full_idx0, int dx_mapped_size1, int al_idx2, int inv_mapping_period3, int *restrict dx_idx_mappings4, long double *restrict der_vec5) {
  return ((al_idx2 < inv_mapping_period3) ? ((dx_idx_mappings4[((inv_mapping_period3 * full_idx0) + al_idx2)] >= 0) ? der_vec5[((dx_mapped_size1 * full_idx0) + dx_idx_mappings4[((inv_mapping_period3 * full_idx0) + al_idx2)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx0, int dx_mapped_size1, int al_idx2, int inv_mapping_period3, int *restrict dx_idx_mappings4, long double *restrict der_vec5) {
  return ((al_idx2 < inv_mapping_period3) ? ((al_idx2 == full_idx0) ? 1.0 : ((dx_idx_mappings4[((inv_mapping_period3 * full_idx0) + al_idx2)] >= 0) ? der_vec5[((dx_mapped_size1 * full_idx0) + dx_idx_mappings4[((inv_mapping_period3 * full_idx0) + al_idx2)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx7, int inv_mapping_period8, int *restrict dx_idx_mappings9, long double *restrict der_vec10) {
  return ((((al_idx7 < inv_mapping_period8) ? dx_idx_mappings9[al_idx7] : -1) >= 0) ? der_vec10[((al_idx7 < inv_mapping_period8) ? dx_idx_mappings9[al_idx7] : -1)] : 0.0);
}

const long double central_grav11[5] = { 0.0e0, 0.0e0, -0.0015242828189020565498e0, 0.0e0, 1.11333333333333329e-05 };

int neptune_rotation_matrix(long double *restrict neptune_rotation_matrix12, long double t13) {
  long double T14;
  T14 = (t13 / 36525.0e0);
  long double N15;
  N15 = (((357.85000000000002274e0 + (52.3160000000000025e0 * T14)) * 3.141592653589793116e0) / 180.0e0);
  long double alpha_016;
  alpha_016 = (((299.36000000000001364e0 + (0.6999999999999999556e0 * sinl(N15))) * 3.141592653589793116e0) / 180.0e0);
  long double delta_017;
  delta_017 = (((43.460000000000000853e0 - (0.5100000000000000089e0 * cosl(N15))) * 3.141592653589793116e0) / 180.0e0);
  long double W18;
  W18 = ((((249.97800000000000864e0 + (541.13977569999997286e0 * t13)) - (0.47999999999999998224e0 * sinl(N15))) * 3.141592653589793116e0) / 180.0e0);
  neptune_rotation_matrix12[0] = ((-sinl(alpha_016) * cosl(W18)) - ((cosl(alpha_016) * sinl(delta_017)) * sinl(W18)));
  neptune_rotation_matrix12[1] = ((sinl(alpha_016) * sinl(W18)) - ((cosl(alpha_016) * sinl(delta_017)) * cosl(W18)));
  neptune_rotation_matrix12[2] = (cosl(alpha_016) * cosl(delta_017));
  neptune_rotation_matrix12[3] = ((cosl(alpha_016) * cosl(W18)) - ((sinl(alpha_016) * sinl(delta_017)) * sinl(W18)));
  neptune_rotation_matrix12[4] = ((-cosl(alpha_016) * sinl(W18)) - ((sinl(alpha_016) * sinl(delta_017)) * cosl(W18)));
  neptune_rotation_matrix12[5] = (cosl(delta_017) * sinl(alpha_016));
  neptune_rotation_matrix12[6] = (cosl(delta_017) * sinl(W18));
  neptune_rotation_matrix12[7] = (cosl(delta_017) * cosl(W18));
  neptune_rotation_matrix12[8] = sinl(delta_017);
  return 0;
}

int nepsatsystem(long double *restrict nepsatsystem19, long double t20, long double *restrict state_and_derivatives21, long double *restrict central_pos22, long double *restrict perturb_gms23, long double *restrict perturb_pos24, long double *restrict sat_gms25) {
  static int ddpotential_dy5697dinitial_mpg26[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddpotential_dy5697dinitial_inv_mpg27[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist32747dinitial_mpg28[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist32747dinitial_inv_mpg29[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dy4126dinitial_mpg30[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dy4126dinitial_inv_mpg31[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_y4276dinitial_mpg32[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_y4276dinitial_inv_mpg33[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r4426dinitial_mpg34[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r4426dinitial_inv_mpg35[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dz3965dinitial_mpg36[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dz3965dinitial_inv_mpg37[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dacc3871dinitial_mpg38[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dacc3871dinitial_inv_mpg39[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dgrav_acc_local3827dinitial_mpg40[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dgrav_acc_local3827dinitial_inv_mpg41[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr2981dinitial_mpg42[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr2981dinitial_inv_mpg43[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r24383dinitial_mpg44[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r24383dinitial_inv_mpg45[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dx4077dinitial_mpg46[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dx4077dinitial_inv_mpg47[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r54499dinitial_mpg48[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r54499dinitial_inv_mpg49[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dy3943dinitial_mpg50[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dy3943dinitial_inv_mpg51[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dstate2201dinitial_mpg52[144] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2201dinitial_inv_mpg53[144] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dz4175dinitial_mpg54[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dz4175dinitial_inv_mpg55[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddpotential_dz5724dinitial_mpg56[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddpotential_dz5724dinitial_inv_mpg57[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dcentral_acc2123dinitial_mpg58[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dcentral_acc2123dinitial_inv_mpg59[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r44474dinitial_mpg60[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r44474dinitial_inv_mpg61[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r34450dinitial_mpg62[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_r34450dinitial_inv_mpg63[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_x4225dinitial_mpg64[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_x4225dinitial_inv_mpg65[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_z4327dinitial_mpg66[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int d_z4327dinitial_inv_mpg67[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dgrav_acc3852dinitial_mpg68[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dgrav_acc3852dinitial_inv_mpg69[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddpotential_dx5670dinitial_mpg70[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddpotential_dx5670dinitial_inv_mpg71[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist23051dinitial_mpg72[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist23051dinitial_inv_mpg73[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist22690dinitial_mpg74[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist22690dinitial_inv_mpg75[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dx3921dinitial_mpg76[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddv_0_dx3921dinitial_inv_mpg77[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dsat_acc2105dinitial_mpg78[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dsat_acc2105dinitial_inv_mpg79[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr3497dinitial_mpg80[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr3497dinitial_inv_mpg81[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dv_03903dinitial_mpg82[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_03903dinitial_inv_mpg83[60] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int ddist23601dinitial_mpg84[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int ddist23601dinitial_inv_mpg85[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr2634dinitial_mpg86[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static int dr2634dinitial_inv_mpg87[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  long double sat_acc89[6] = { 0.0 };
  long double dsat_acc2105dinitial_der88[72] = { 0.0 };
  long double central_acc91[3] = { 0.0 };
  long double dcentral_acc2123dinitial_der90[36] = { 0.0 };
  long double state_derivatives_initial92[144] = { 0.0 };
  long double state94[12] = { 0.0 };
  long double dstate2201dinitial_der93[144] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 144; slice_idx++) {
      state_derivatives_initial92[(0 + slice_idx)] = state_and_derivatives21[(12 + slice_idx)];
    }
  }
  for (int slice_idx = 0; slice_idx < 12; slice_idx++) {
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
      if ((mappings_full_idx_symbol >= 144)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dstate2201dinitial_mpg52[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dstate2201dinitial_der93[((12 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 12; slice_idx++) {
    state94[(0 + slice_idx)] = state_and_derivatives21[(0 + slice_idx)];
  }
  for (int i95 = 0; i95 < 12; i95++) {
    for (int j96 = 0; j96 < 12; j96++) {
      {
        int mapped_idx;
        mapped_idx = ((j96 < 12) ? dstate2201dinitial_inv_mpg53[((i95 * 12) + j96)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2201dinitial_der93[((i95 * 12) + mapped_idx)] = state_derivatives_initial92[(((i95 * 2) * 6) + j96)];
        }
      }
    }
  }
  for (int i97 = 0; i97 < 2; i97++) {
    long double r99[3] = { 0.0 };
    long double dr2634dinitial_der98[36] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr2634dinitial_mpg86[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dr2634dinitial_der98[((12 * (0 + slice_idx)) + mapped_idx)] = get_dfdx_cell(((i97 * 6) + slice_idx), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r99[(0 + slice_idx)] = state94[((i97 * 6) + slice_idx)];
    }
    long double dist2101;
    long double ddist22690dinitial_der100[12] = { 0.0 };
    for (int loop_var102 = 0; loop_var102 < 12; loop_var102++) {
      int mapped_idx;
      mapped_idx = loop_var102;
      int al_index_name_symbol;
      al_index_name_symbol = ddist22690dinitial_mpg74[loop_var102];
      ddist22690dinitial_der100[mapped_idx] = ((((r99[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98)) + (r99[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98))) + ((r99[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98)) + (r99[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98)))) + ((r99[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98)) + (r99[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98))));
    }
    dist2101 = (((r99[0] * r99[0]) + (r99[1] * r99[1])) + (r99[2] * r99[2]));
    long double dist3104;
    long double ddist32747dinitial_der103[12] = { 0.0 };
    for (int loop_var105 = 0; loop_var105 < 12; loop_var105++) {
      int mapped_idx;
      mapped_idx = loop_var105;
      int al_index_name_symbol;
      al_index_name_symbol = ddist32747dinitial_mpg28[loop_var105];
      ddist32747dinitial_der103[mapped_idx] = ((dist2101 * (0.5 * (sqrtl((1.0 / dist2101)) * get_dfdx_var(al_index_name_symbol, 12, ddist22690dinitial_inv_mpg75, ddist22690dinitial_der100)))) + (sqrtl(dist2101) * get_dfdx_var(al_index_name_symbol, 12, ddist22690dinitial_inv_mpg75, ddist22690dinitial_der100)));
    }
    dist3104 = (dist2101 * sqrtl(dist2101));
    for (int slice_idx = 0; slice_idx < (((i97 * 3) + 3) - (i97 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * ((i97 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2105dinitial_mpg78[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dsat_acc2105dinitial_der88[((12 * ((i97 * 3) + slice_idx)) + mapped_idx)] = ((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98) * -1.5240390322546755564e-08) * dist3104) - (get_dfdx_var(al_index_name_symbol, 12, ddist32747dinitial_inv_mpg29, ddist32747dinitial_der103) * (r99[(0 + slice_idx)] * -1.5240390322546755564e-08))) / (dist3104 * dist3104));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i97 * 3) + 3) - (i97 * 3)); slice_idx++) {
      sat_acc89[((i97 * 3) + slice_idx)] = ((-1.5240390322546755564e-08 * r99[(0 + slice_idx)]) / dist3104);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2123dinitial_mpg58[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dcentral_acc2123dinitial_der90[((12 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dcentral_acc2123dinitial_inv_mpg59, dcentral_acc2123dinitial_der90) + ((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr2634dinitial_inv_mpg87, dr2634dinitial_der98) * sat_gms25[i97]) * dist3104) - (get_dfdx_var(al_index_name_symbol, 12, ddist32747dinitial_inv_mpg29, ddist32747dinitial_der103) * (r99[(0 + slice_idx)] * sat_gms25[i97]))) / (dist3104 * dist3104)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc91[(0 + slice_idx)] = (central_acc91[(0 + slice_idx)] + ((sat_gms25[i97] * r99[(0 + slice_idx)]) / dist3104));
    }
  }
  for (int i106 = 0; i106 < 4; i106++) {
    long double r108[3] = { 0.0 };
    long double dr2981dinitial_der107[36] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr2981dinitial_mpg42[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dr2981dinitial_der107[((12 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r108[(0 + slice_idx)] = (perturb_pos24[((i106 * 3) + slice_idx)] - central_pos22[(0 + slice_idx)]);
    }
    long double dist2110;
    long double ddist23051dinitial_der109[12] = { 0.0 };
    for (int loop_var111 = 0; loop_var111 < 12; loop_var111++) {
      int mapped_idx;
      mapped_idx = loop_var111;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23051dinitial_mpg72[loop_var111];
      ddist23051dinitial_der109[mapped_idx] = ((((r108[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107))) + ((r108[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)))) + ((r108[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107))));
    }
    dist2110 = (((r108[0] * r108[0]) + (r108[1] * r108[1])) + (r108[2] * r108[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2123dinitial_mpg58[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dcentral_acc2123dinitial_der90[((12 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dcentral_acc2123dinitial_inv_mpg59, dcentral_acc2123dinitial_der90) + (((((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107) * perturb_gms23[i106]) * dist2110) - (get_dfdx_var(al_index_name_symbol, 12, ddist23051dinitial_inv_mpg73, ddist23051dinitial_der109) * (r108[(0 + slice_idx)] * perturb_gms23[i106]))) / (dist2110 * dist2110)) * sqrtl(dist2110)) - ((0.5 * (sqrtl((1.0 / dist2110)) * get_dfdx_var(al_index_name_symbol, 12, ddist23051dinitial_inv_mpg73, ddist23051dinitial_der109))) * ((r108[(0 + slice_idx)] * perturb_gms23[i106]) / dist2110))) / (sqrtl(dist2110) * sqrtl(dist2110))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc91[(0 + slice_idx)] = (central_acc91[(0 + slice_idx)] + (((perturb_gms23[i106] * r108[(0 + slice_idx)]) / dist2110) / sqrtl(dist2110)));
    }
    for (int j112 = 0; j112 < 2; j112++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 36)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr2981dinitial_mpg42[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
              break;
            } else {
              dr2981dinitial_der107[((12 * (0 + slice_idx)) + mapped_idx)] = -get_dfdx_cell(((j112 * 6) + slice_idx), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r108[(0 + slice_idx)] = (perturb_pos24[((i106 * 3) + slice_idx)] - (state94[((j112 * 6) + slice_idx)] + central_pos22[(0 + slice_idx)]));
      }
      for (int loop_var113 = 0; loop_var113 < 12; loop_var113++) {
        int mapped_idx;
        mapped_idx = loop_var113;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23051dinitial_mpg72[loop_var113];
        ddist23051dinitial_der109[mapped_idx] = ((((r108[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107))) + ((r108[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)))) + ((r108[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107)) + (r108[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107))));
      }
      dist2110 = (((r108[0] * r108[0]) + (r108[1] * r108[1])) + (r108[2] * r108[2]));
      for (int slice_idx = 0; slice_idx < (((j112 * 3) + 3) - (j112 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (12 * ((j112 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2105dinitial_mpg78[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
              break;
            } else {
              dsat_acc2105dinitial_der88[((12 * ((j112 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j112 * 3) + slice_idx), 12, al_index_name_symbol, 12, dsat_acc2105dinitial_inv_mpg79, dsat_acc2105dinitial_der88) + (((((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr2981dinitial_inv_mpg43, dr2981dinitial_der107) * perturb_gms23[i106]) * dist2110) - (get_dfdx_var(al_index_name_symbol, 12, ddist23051dinitial_inv_mpg73, ddist23051dinitial_der109) * (r108[(0 + slice_idx)] * perturb_gms23[i106]))) / (dist2110 * dist2110)) * sqrtl(dist2110)) - ((0.5 * (sqrtl((1.0 / dist2110)) * get_dfdx_var(al_index_name_symbol, 12, ddist23051dinitial_inv_mpg73, ddist23051dinitial_der109))) * ((r108[(0 + slice_idx)] * perturb_gms23[i106]) / dist2110))) / (sqrtl(dist2110) * sqrtl(dist2110))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j112 * 3) + 3) - (j112 * 3)); slice_idx++) {
        sat_acc89[((j112 * 3) + slice_idx)] = (sat_acc89[((j112 * 3) + slice_idx)] + (((perturb_gms23[i106] * r108[(0 + slice_idx)]) / dist2110) / sqrtl(dist2110)));
      }
    }
  }
  for (int i114 = 1; i114 < 2; i114++) {
    long double r116[3] = { 0.0 };
    long double dr3497dinitial_der115[36] = { 0.0 };
    for (int j117 = 0; j117 < i114; j117++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 36)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3497dinitial_mpg80[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
              break;
            } else {
              dr3497dinitial_der115[((12 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j117 * 6) + slice_idx), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93) - get_dfdx_cell(((i114 * 6) + slice_idx), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r116[(0 + slice_idx)] = (state94[((j117 * 6) + slice_idx)] - state94[((i114 * 6) + slice_idx)]);
      }
      long double dist2119;
      long double ddist23601dinitial_der118[12] = { 0.0 };
      for (int loop_var120 = 0; loop_var120 < 12; loop_var120++) {
        int mapped_idx;
        mapped_idx = loop_var120;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23601dinitial_mpg84[loop_var120];
        ddist23601dinitial_der118[mapped_idx] = ((((r116[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115)) + (r116[0] * get_dfdx_cell(0, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115))) + ((r116[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115)) + (r116[1] * get_dfdx_cell(1, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115)))) + ((r116[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115)) + (r116[2] * get_dfdx_cell(2, 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115))));
      }
      dist2119 = (((r116[0] * r116[0]) + (r116[1] * r116[1])) + (r116[2] * r116[2]));
      for (int slice_idx = 0; slice_idx < (((i114 * 3) + 3) - (i114 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (12 * ((i114 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2105dinitial_mpg78[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
              break;
            } else {
              dsat_acc2105dinitial_der88[((12 * ((i114 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i114 * 3) + slice_idx), 12, al_index_name_symbol, 12, dsat_acc2105dinitial_inv_mpg79, dsat_acc2105dinitial_der88) + (((((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115) * sat_gms25[j117]) * dist2119) - (get_dfdx_var(al_index_name_symbol, 12, ddist23601dinitial_inv_mpg85, ddist23601dinitial_der118) * (r116[(0 + slice_idx)] * sat_gms25[j117]))) / (dist2119 * dist2119)) * sqrtl(dist2119)) - ((0.5 * (sqrtl((1.0 / dist2119)) * get_dfdx_var(al_index_name_symbol, 12, ddist23601dinitial_inv_mpg85, ddist23601dinitial_der118))) * ((r116[(0 + slice_idx)] * sat_gms25[j117]) / dist2119))) / (sqrtl(dist2119) * sqrtl(dist2119))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i114 * 3) + 3) - (i114 * 3)); slice_idx++) {
        sat_acc89[((i114 * 3) + slice_idx)] = (sat_acc89[((i114 * 3) + slice_idx)] + (((sat_gms25[j117] * r116[(0 + slice_idx)]) / dist2119) / sqrtl(dist2119)));
      }
      for (int slice_idx = 0; slice_idx < (((j117 * 3) + 3) - (j117 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (12 * ((j117 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2105dinitial_mpg78[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
              break;
            } else {
              dsat_acc2105dinitial_der88[((12 * ((j117 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j117 * 3) + slice_idx), 12, al_index_name_symbol, 12, dsat_acc2105dinitial_inv_mpg79, dsat_acc2105dinitial_der88) - (((((((get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dr3497dinitial_inv_mpg81, dr3497dinitial_der115) * sat_gms25[i114]) * dist2119) - (get_dfdx_var(al_index_name_symbol, 12, ddist23601dinitial_inv_mpg85, ddist23601dinitial_der118) * (r116[(0 + slice_idx)] * sat_gms25[i114]))) / (dist2119 * dist2119)) * sqrtl(dist2119)) - ((0.5 * (sqrtl((1.0 / dist2119)) * get_dfdx_var(al_index_name_symbol, 12, ddist23601dinitial_inv_mpg85, ddist23601dinitial_der118))) * ((r116[(0 + slice_idx)] * sat_gms25[i114]) / dist2119))) / (sqrtl(dist2119) * sqrtl(dist2119))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j117 * 3) + 3) - (j117 * 3)); slice_idx++) {
        sat_acc89[((j117 * 3) + slice_idx)] = (sat_acc89[((j117 * 3) + slice_idx)] - (((sat_gms25[i114] * r116[(0 + slice_idx)]) / dist2119) / sqrtl(dist2119)));
      }
    }
  }
  long double grav_acc_local122[3] = { 0.0 };
  long double dgrav_acc_local3827dinitial_der121[36] = { 0.0 };
  long double grav_acc124[3] = { 0.0 };
  long double dgrav_acc3852dinitial_der123[36] = { 0.0 };
  long double acc126[3] = { 0.0 };
  long double dacc3871dinitial_der125[36] = { 0.0 };
  long double rot127[9] = { 0.0 };
  long double v_0129[5] = { 0.0 };
  long double dv_03903dinitial_der128[60] = { 0.0 };
  long double dv_0_dx131[5] = { 0.0 };
  long double ddv_0_dx3921dinitial_der130[60] = { 0.0 };
  long double dv_0_dy133[5] = { 0.0 };
  long double ddv_0_dy3943dinitial_der132[60] = { 0.0 };
  long double dv_0_dz135[5] = { 0.0 };
  long double ddv_0_dz3965dinitial_der134[60] = { 0.0 };
  {
    long double neptune_rotation_matrix3985137[9] = { 0.0 };
    long double _inl_neptune_rotation_matrix_T139;
    _inl_neptune_rotation_matrix_T139 = (t20 / 36525.0e0);
    long double _inl_neptune_rotation_matrix_N140;
    _inl_neptune_rotation_matrix_N140 = (((357.85000000000002274e0 + (52.3160000000000025e0 * _inl_neptune_rotation_matrix_T139)) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_alpha_0141;
    _inl_neptune_rotation_matrix_alpha_0141 = (((299.36000000000001364e0 + (0.6999999999999999556e0 * sinl(_inl_neptune_rotation_matrix_N140))) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_delta_0142;
    _inl_neptune_rotation_matrix_delta_0142 = (((43.460000000000000853e0 - (0.5100000000000000089e0 * cosl(_inl_neptune_rotation_matrix_N140))) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_W143;
    _inl_neptune_rotation_matrix_W143 = ((((249.97800000000000864e0 + (541.13977569999997286e0 * t20)) - (0.47999999999999998224e0 * sinl(_inl_neptune_rotation_matrix_N140))) * 3.141592653589793116e0) / 180.0e0);
    neptune_rotation_matrix3985137[0] = ((-sinl(_inl_neptune_rotation_matrix_alpha_0141) * cosl(_inl_neptune_rotation_matrix_W143)) - ((cosl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_delta_0142)) * sinl(_inl_neptune_rotation_matrix_W143)));
    neptune_rotation_matrix3985137[1] = ((sinl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_W143)) - ((cosl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_delta_0142)) * cosl(_inl_neptune_rotation_matrix_W143)));
    neptune_rotation_matrix3985137[2] = (cosl(_inl_neptune_rotation_matrix_alpha_0141) * cosl(_inl_neptune_rotation_matrix_delta_0142));
    neptune_rotation_matrix3985137[3] = ((cosl(_inl_neptune_rotation_matrix_alpha_0141) * cosl(_inl_neptune_rotation_matrix_W143)) - ((sinl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_delta_0142)) * sinl(_inl_neptune_rotation_matrix_W143)));
    neptune_rotation_matrix3985137[4] = ((-cosl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_W143)) - ((sinl(_inl_neptune_rotation_matrix_alpha_0141) * sinl(_inl_neptune_rotation_matrix_delta_0142)) * cosl(_inl_neptune_rotation_matrix_W143)));
    neptune_rotation_matrix3985137[5] = (cosl(_inl_neptune_rotation_matrix_delta_0142) * sinl(_inl_neptune_rotation_matrix_alpha_0141));
    neptune_rotation_matrix3985137[6] = (cosl(_inl_neptune_rotation_matrix_delta_0142) * sinl(_inl_neptune_rotation_matrix_W143));
    neptune_rotation_matrix3985137[7] = (cosl(_inl_neptune_rotation_matrix_delta_0142) * cosl(_inl_neptune_rotation_matrix_W143));
    neptune_rotation_matrix3985137[8] = sinl(_inl_neptune_rotation_matrix_delta_0142);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot127[(0 + slice_idx)] = neptune_rotation_matrix3985137[(0 + slice_idx)];
    }
  }
  for (int i144 = 0; i144 < 2; i144++) {
    long double x146;
    long double dx4077dinitial_der145[12] = { 0.0 };
    for (int loop_var147 = 0; loop_var147 < 12; loop_var147++) {
      int mapped_idx;
      mapped_idx = loop_var147;
      int al_index_name_symbol;
      al_index_name_symbol = dx4077dinitial_mpg46[loop_var147];
      dx4077dinitial_der145[mapped_idx] = (get_dfdx_cell(((i144 * 6) + 0), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93) / 0.00016861871015922155108e0);
    }
    x146 = (state94[((i144 * 6) + 0)] / 0.00016861871015922155108e0);
    long double y149;
    long double dy4126dinitial_der148[12] = { 0.0 };
    for (int loop_var150 = 0; loop_var150 < 12; loop_var150++) {
      int mapped_idx;
      mapped_idx = loop_var150;
      int al_index_name_symbol;
      al_index_name_symbol = dy4126dinitial_mpg30[loop_var150];
      dy4126dinitial_der148[mapped_idx] = (get_dfdx_cell(((i144 * 6) + 1), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93) / 0.00016861871015922155108e0);
    }
    y149 = (state94[((i144 * 6) + 1)] / 0.00016861871015922155108e0);
    long double z152;
    long double dz4175dinitial_der151[12] = { 0.0 };
    for (int loop_var153 = 0; loop_var153 < 12; loop_var153++) {
      int mapped_idx;
      mapped_idx = loop_var153;
      int al_index_name_symbol;
      al_index_name_symbol = dz4175dinitial_mpg54[loop_var153];
      dz4175dinitial_der151[mapped_idx] = (get_dfdx_cell(((i144 * 6) + 2), 12, al_index_name_symbol, 12, dstate2201dinitial_inv_mpg53, dstate2201dinitial_der93) / 0.00016861871015922155108e0);
    }
    z152 = (state94[((i144 * 6) + 2)] / 0.00016861871015922155108e0);
    long double _x155;
    long double d_x4225dinitial_der154[12] = { 0.0 };
    for (int loop_var156 = 0; loop_var156 < 12; loop_var156++) {
      int mapped_idx;
      mapped_idx = loop_var156;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4225dinitial_mpg64[loop_var156];
      d_x4225dinitial_der154[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 12, dx4077dinitial_inv_mpg47, dx4077dinitial_der145) * rot127[0]) + (get_dfdx_var(al_index_name_symbol, 12, dy4126dinitial_inv_mpg31, dy4126dinitial_der148) * rot127[3])) + (get_dfdx_var(al_index_name_symbol, 12, dz4175dinitial_inv_mpg55, dz4175dinitial_der151) * rot127[6]));
    }
    _x155 = (((rot127[0] * x146) + (rot127[3] * y149)) + (rot127[6] * z152));
    long double _y158;
    long double d_y4276dinitial_der157[12] = { 0.0 };
    for (int loop_var159 = 0; loop_var159 < 12; loop_var159++) {
      int mapped_idx;
      mapped_idx = loop_var159;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4276dinitial_mpg32[loop_var159];
      d_y4276dinitial_der157[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 12, dx4077dinitial_inv_mpg47, dx4077dinitial_der145) * rot127[1]) + (get_dfdx_var(al_index_name_symbol, 12, dy4126dinitial_inv_mpg31, dy4126dinitial_der148) * rot127[4])) + (get_dfdx_var(al_index_name_symbol, 12, dz4175dinitial_inv_mpg55, dz4175dinitial_der151) * rot127[7]));
    }
    _y158 = (((rot127[1] * x146) + (rot127[4] * y149)) + (rot127[7] * z152));
    long double _z161;
    long double d_z4327dinitial_der160[12] = { 0.0 };
    for (int loop_var162 = 0; loop_var162 < 12; loop_var162++) {
      int mapped_idx;
      mapped_idx = loop_var162;
      int al_index_name_symbol;
      al_index_name_symbol = d_z4327dinitial_mpg66[loop_var162];
      d_z4327dinitial_der160[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 12, dx4077dinitial_inv_mpg47, dx4077dinitial_der145) * rot127[2]) + (get_dfdx_var(al_index_name_symbol, 12, dy4126dinitial_inv_mpg31, dy4126dinitial_der148) * rot127[5])) + (get_dfdx_var(al_index_name_symbol, 12, dz4175dinitial_inv_mpg55, dz4175dinitial_der151) * rot127[8]));
    }
    _z161 = (((rot127[2] * x146) + (rot127[5] * y149)) + (rot127[8] * z152));
    long double _r2164;
    long double d_r24383dinitial_der163[12] = { 0.0 };
    for (int loop_var165 = 0; loop_var165 < 12; loop_var165++) {
      int mapped_idx;
      mapped_idx = loop_var165;
      int al_index_name_symbol;
      al_index_name_symbol = d_r24383dinitial_mpg44[loop_var165];
      d_r24383dinitial_der163[mapped_idx] = ((((_x155 * get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154)) + (_x155 * get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154))) + ((_y158 * get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157)) + (_y158 * get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157)))) + ((_z161 * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160))));
    }
    _r2164 = (((_x155 * _x155) + (_y158 * _y158)) + (_z161 * _z161));
    long double _r167;
    long double d_r4426dinitial_der166[12] = { 0.0 };
    for (int loop_var168 = 0; loop_var168 < 12; loop_var168++) {
      int mapped_idx;
      mapped_idx = loop_var168;
      int al_index_name_symbol;
      al_index_name_symbol = d_r4426dinitial_mpg34[loop_var168];
      d_r4426dinitial_der166[mapped_idx] = (0.5 * (sqrtl((1.0 / _r2164)) * get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163)));
    }
    _r167 = sqrtl(_r2164);
    long double _r3170;
    long double d_r34450dinitial_der169[12] = { 0.0 };
    for (int loop_var171 = 0; loop_var171 < 12; loop_var171++) {
      int mapped_idx;
      mapped_idx = loop_var171;
      int al_index_name_symbol;
      al_index_name_symbol = d_r34450dinitial_mpg62[loop_var171];
      d_r34450dinitial_der169[mapped_idx] = ((_r167 * get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163)) + (_r2164 * get_dfdx_var(al_index_name_symbol, 12, d_r4426dinitial_inv_mpg35, d_r4426dinitial_der166)));
    }
    _r3170 = (_r167 * _r2164);
    long double _r4173;
    long double d_r44474dinitial_der172[12] = { 0.0 };
    for (int loop_var174 = 0; loop_var174 < 12; loop_var174++) {
      int mapped_idx;
      mapped_idx = loop_var174;
      int al_index_name_symbol;
      al_index_name_symbol = d_r44474dinitial_mpg60[loop_var174];
      d_r44474dinitial_der172[mapped_idx] = ((_r2164 * get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163)) + (_r2164 * get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163)));
    }
    _r4173 = (_r2164 * _r2164);
    long double _r5176;
    long double d_r54499dinitial_der175[12] = { 0.0 };
    for (int loop_var177 = 0; loop_var177 < 12; loop_var177++) {
      int mapped_idx;
      mapped_idx = loop_var177;
      int al_index_name_symbol;
      al_index_name_symbol = d_r54499dinitial_mpg48[loop_var177];
      d_r54499dinitial_der175[mapped_idx] = ((_r4173 * get_dfdx_var(al_index_name_symbol, 12, d_r4426dinitial_inv_mpg35, d_r4426dinitial_der166)) + (_r167 * get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172)));
    }
    _r5176 = (_r4173 * _r167);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 0));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_03903dinitial_mpg82[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dv_03903dinitial_der128[((12 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r167 * _r167)) * get_dfdx_var(al_index_name_symbol, 12, d_r4426dinitial_inv_mpg35, d_r4426dinitial_der166))));
        }
      }
    }
    v_0129[0] = (1.0e0 / _r167);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 1));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_03903dinitial_mpg82[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dv_03903dinitial_der128[((12 * 1) + mapped_idx)] = ((((get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * 1.7320508075688772936e0) * _r3170) - (get_dfdx_var(al_index_name_symbol, 12, d_r34450dinitial_inv_mpg63, d_r34450dinitial_der169) * (_z161 * 1.7320508075688772936e0))) / (_r3170 * _r3170));
        }
      }
    }
    v_0129[1] = ((1.7320508075688772936e0 * _z161) / _r3170);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 0));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx3921dinitial_mpg76[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dx3921dinitial_der130[((12 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154) * _r3170) - (get_dfdx_var(al_index_name_symbol, 12, d_r34450dinitial_inv_mpg63, d_r34450dinitial_der169) * -_x155)) / (_r3170 * _r3170));
        }
      }
    }
    dv_0_dx131[0] = (-_x155 / _r3170);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 0));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy3943dinitial_mpg50[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dy3943dinitial_der132[((12 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157) * _r3170) - (get_dfdx_var(al_index_name_symbol, 12, d_r34450dinitial_inv_mpg63, d_r34450dinitial_der169) * -_y158)) / (_r3170 * _r3170));
        }
      }
    }
    dv_0_dy133[0] = (-_y158 / _r3170);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 0));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz3965dinitial_mpg36[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dz3965dinitial_der134[((12 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * _r3170) - (get_dfdx_var(al_index_name_symbol, 12, d_r34450dinitial_inv_mpg63, d_r34450dinitial_der169) * -_z161)) / (_r3170 * _r3170));
        }
      }
    }
    dv_0_dz135[0] = (-_z161 / _r3170);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 1));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx3921dinitial_mpg76[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dx3921dinitial_der130[((12 * 1) + mapped_idx)] = ((((((_z161 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154)) + (_x155 * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * -5.1961524227066318805e0))) * _r5176) - (get_dfdx_var(al_index_name_symbol, 12, d_r54499dinitial_inv_mpg49, d_r54499dinitial_der175) * ((_z161 * -5.1961524227066318805e0) * _x155))) / (_r5176 * _r5176));
        }
      }
    }
    dv_0_dx131[1] = (((-5.1961524227066318805e0 * _z161) * _x155) / _r5176);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 1));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy3943dinitial_mpg50[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dy3943dinitial_der132[((12 * 1) + mapped_idx)] = ((((((_z161 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157)) + (_y158 * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * -5.1961524227066318805e0))) * _r5176) - (get_dfdx_var(al_index_name_symbol, 12, d_r54499dinitial_inv_mpg49, d_r54499dinitial_der175) * ((_z161 * -5.1961524227066318805e0) * _y158))) / (_r5176 * _r5176));
        }
      }
    }
    dv_0_dy133[1] = (((-5.1961524227066318805e0 * _z161) * _y158) / _r5176);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 1));
      if ((mappings_full_idx_symbol >= 60)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz3965dinitial_mpg36[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          ddv_0_dz3965dinitial_der134[((12 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r3170 * _r3170)) * get_dfdx_var(al_index_name_symbol, 12, d_r34450dinitial_inv_mpg63, d_r34450dinitial_der169)))) - ((((((_z161 * 3.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * 3.0e0))) * _r5176) - (get_dfdx_var(al_index_name_symbol, 12, d_r54499dinitial_inv_mpg49, d_r54499dinitial_der175) * ((_z161 * 3.0e0) * _z161))) / (_r5176 * _r5176))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz135[1] = (1.7320508075688772936e0 * ((1.0e0 / _r3170) - (((3.0e0 * _z161) * _z161) / _r5176)));
    for (int n178 = 2; n178 < 5; n178++) {
      long double coef1179;
      coef1179 = sqrtl(((((2.0e0 * ((long double) n178)) - 1.0e0) * ((long double) ((2 * n178) + 1))) / ((long double) (n178 * n178))));
      long double coef2180;
      coef2180 = sqrtl(((((((long double) n178) - 1.0e0) * ((long double) (n178 - 1))) * ((long double) ((2 * n178) + 1))) / ((long double) ((n178 * n178) * ((2 * n178) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * n178));
        if ((mappings_full_idx_symbol >= 60)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_03903dinitial_mpg82[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dv_03903dinitial_der128[((12 * n178) + mapped_idx)] = (((((((v_0129[(n178 - 1)] * coef1179) * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * (get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * coef1179))) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * ((v_0129[(n178 - 1)] * coef1179) * _z161))) / (_r2164 * _r2164)) - ((((get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * coef2180) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * (v_0129[(n178 - 2)] * coef2180))) / (_r2164 * _r2164)));
          }
        }
      }
      v_0129[n178] = ((((coef1179 * v_0129[(n178 - 1)]) * _z161) / _r2164) - ((coef2180 * v_0129[(n178 - 2)]) / _r2164));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * n178));
        if ((mappings_full_idx_symbol >= 60)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx3921dinitial_mpg76[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            ddv_0_dx3921dinitial_der130[((12 * n178) + mapped_idx)] = ((((_z161 * coef1179) * ((((get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, ddv_0_dx3921dinitial_inv_mpg77, ddv_0_dx3921dinitial_der130) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * dv_0_dx131[(n178 - 1)])) / (_r2164 * _r2164)) - ((((((v_0129[(n178 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154)) + (_x155 * (get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((v_0129[(n178 - 1)] * 2.0e0) * _x155))) / (_r4173 * _r4173)))) + (((dv_0_dx131[(n178 - 1)] / _r2164) - (((v_0129[(n178 - 1)] * 2.0e0) * _x155) / _r4173)) * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * coef1179))) - (((((get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, ddv_0_dx3921dinitial_inv_mpg77, ddv_0_dx3921dinitial_der130) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * dv_0_dx131[(n178 - 2)])) / (_r2164 * _r2164)) - ((((((v_0129[(n178 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_x4225dinitial_inv_mpg65, d_x4225dinitial_der154)) + (_x155 * (get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((v_0129[(n178 - 2)] * 2.0e0) * _x155))) / (_r4173 * _r4173))) * coef2180));
          }
        }
      }
      dv_0_dx131[n178] = (((coef1179 * _z161) * ((dv_0_dx131[(n178 - 1)] / _r2164) - (((v_0129[(n178 - 1)] * 2.0e0) * _x155) / _r4173))) - (coef2180 * ((dv_0_dx131[(n178 - 2)] / _r2164) - (((v_0129[(n178 - 2)] * 2.0e0) * _x155) / _r4173))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * n178));
        if ((mappings_full_idx_symbol >= 60)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy3943dinitial_mpg50[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            ddv_0_dy3943dinitial_der132[((12 * n178) + mapped_idx)] = ((((_z161 * coef1179) * ((((get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, ddv_0_dy3943dinitial_inv_mpg51, ddv_0_dy3943dinitial_der132) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * dv_0_dy133[(n178 - 1)])) / (_r2164 * _r2164)) - ((((((v_0129[(n178 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157)) + (_y158 * (get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((v_0129[(n178 - 1)] * 2.0e0) * _y158))) / (_r4173 * _r4173)))) + (((dv_0_dy133[(n178 - 1)] / _r2164) - (((v_0129[(n178 - 1)] * 2.0e0) * _y158) / _r4173)) * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * coef1179))) - (((((get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, ddv_0_dy3943dinitial_inv_mpg51, ddv_0_dy3943dinitial_der132) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * dv_0_dy133[(n178 - 2)])) / (_r2164 * _r2164)) - ((((((v_0129[(n178 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_y4276dinitial_inv_mpg33, d_y4276dinitial_der157)) + (_y158 * (get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((v_0129[(n178 - 2)] * 2.0e0) * _y158))) / (_r4173 * _r4173))) * coef2180));
          }
        }
      }
      dv_0_dy133[n178] = (((coef1179 * _z161) * ((dv_0_dy133[(n178 - 1)] / _r2164) - (((v_0129[(n178 - 1)] * 2.0e0) * _y158) / _r4173))) - (coef2180 * ((dv_0_dy133[(n178 - 2)] / _r2164) - (((v_0129[(n178 - 2)] * 2.0e0) * _y158) / _r4173))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * n178));
        if ((mappings_full_idx_symbol >= 60)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz3965dinitial_mpg36[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            ddv_0_dz3965dinitial_der134[((12 * n178) + mapped_idx)] = ((((((((dv_0_dz135[(n178 - 1)] * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, ddv_0_dz3965dinitial_inv_mpg37, ddv_0_dz3965dinitial_der134))) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * (dv_0_dz135[(n178 - 1)] * _z161))) / (_r2164 * _r2164)) + ((v_0129[(n178 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r2164 * _r2164)) * get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163)))) - ((((((_z161 * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * (get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((_z161 * 2.0e0) * _z161))) / (_r4173 * _r4173)))) + (((1.0e0 / _r2164) - (((_z161 * 2.0e0) * _z161) / _r4173)) * get_dfdx_cell((n178 - 1), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128)))) * coef1179) - (((((get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, ddv_0_dz3965dinitial_inv_mpg37, ddv_0_dz3965dinitial_der134) * _r2164) - (get_dfdx_var(al_index_name_symbol, 12, d_r24383dinitial_inv_mpg45, d_r24383dinitial_der163) * dv_0_dz135[(n178 - 2)])) / (_r2164 * _r2164)) - ((((((v_0129[(n178 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 12, d_z4327dinitial_inv_mpg67, d_z4327dinitial_der160)) + (_z161 * (get_dfdx_cell((n178 - 2), 12, al_index_name_symbol, 12, dv_03903dinitial_inv_mpg83, dv_03903dinitial_der128) * 2.0e0))) * _r4173) - (get_dfdx_var(al_index_name_symbol, 12, d_r44474dinitial_inv_mpg61, d_r44474dinitial_der172) * ((v_0129[(n178 - 2)] * 2.0e0) * _z161))) / (_r4173 * _r4173))) * coef2180));
          }
        }
      }
      dv_0_dz135[n178] = ((coef1179 * (((dv_0_dz135[(n178 - 1)] * _z161) / _r2164) + (v_0129[(n178 - 1)] * ((1.0e0 / _r2164) - (((2.0e0 * _z161) * _z161) / _r4173))))) - (coef2180 * ((dv_0_dz135[(n178 - 2)] / _r2164) - (((v_0129[(n178 - 2)] * 2.0e0) * _z161) / _r4173))));
    }
    long double dpotential_dx182;
    long double ddpotential_dx5670dinitial_der181[12] = { 0.0 };
    for (int loop_var183 = 0; loop_var183 < 12; loop_var183++) {
      int mapped_idx;
      mapped_idx = loop_var183;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx5670dinitial_mpg70[loop_var183];
      ddpotential_dx5670dinitial_der181[mapped_idx] = 0.0;
    }
    dpotential_dx182 = ((long double) 0);
    long double dpotential_dy185;
    long double ddpotential_dy5697dinitial_der184[12] = { 0.0 };
    for (int loop_var186 = 0; loop_var186 < 12; loop_var186++) {
      int mapped_idx;
      mapped_idx = loop_var186;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy5697dinitial_mpg26[loop_var186];
      ddpotential_dy5697dinitial_der184[mapped_idx] = 0.0;
    }
    dpotential_dy185 = ((long double) 0);
    long double dpotential_dz188;
    long double ddpotential_dz5724dinitial_der187[12] = { 0.0 };
    for (int loop_var189 = 0; loop_var189 < 12; loop_var189++) {
      int mapped_idx;
      mapped_idx = loop_var189;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz5724dinitial_mpg56[loop_var189];
      ddpotential_dz5724dinitial_der187[mapped_idx] = 0.0;
    }
    dpotential_dz188 = ((long double) 0);
    for (int n190 = 2; n190 < 5; n190++) {
      for (int loop_var191 = 0; loop_var191 < 12; loop_var191++) {
        int mapped_idx;
        mapped_idx = loop_var191;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx5670dinitial_mpg70[loop_var191];
        ddpotential_dx5670dinitial_der181[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 12, ddpotential_dx5670dinitial_inv_mpg71, ddpotential_dx5670dinitial_der181) + (get_dfdx_cell(n190, 12, al_index_name_symbol, 12, ddv_0_dx3921dinitial_inv_mpg77, ddv_0_dx3921dinitial_der130) * central_grav11[n190]));
      }
      dpotential_dx182 = (dpotential_dx182 + (dv_0_dx131[n190] * central_grav11[n190]));
      for (int loop_var192 = 0; loop_var192 < 12; loop_var192++) {
        int mapped_idx;
        mapped_idx = loop_var192;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy5697dinitial_mpg26[loop_var192];
        ddpotential_dy5697dinitial_der184[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 12, ddpotential_dy5697dinitial_inv_mpg27, ddpotential_dy5697dinitial_der184) + (get_dfdx_cell(n190, 12, al_index_name_symbol, 12, ddv_0_dy3943dinitial_inv_mpg51, ddv_0_dy3943dinitial_der132) * central_grav11[n190]));
      }
      dpotential_dy185 = (dpotential_dy185 + (dv_0_dy133[n190] * central_grav11[n190]));
      for (int loop_var193 = 0; loop_var193 < 12; loop_var193++) {
        int mapped_idx;
        mapped_idx = loop_var193;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz5724dinitial_mpg56[loop_var193];
        ddpotential_dz5724dinitial_der187[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 12, ddpotential_dz5724dinitial_inv_mpg57, ddpotential_dz5724dinitial_der187) + (get_dfdx_cell(n190, 12, al_index_name_symbol, 12, ddv_0_dz3965dinitial_inv_mpg37, ddv_0_dz3965dinitial_der134) * central_grav11[n190]));
      }
      dpotential_dz188 = (dpotential_dz188 + (dv_0_dz135[n190] * central_grav11[n190]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 0));
      if ((mappings_full_idx_symbol >= 36)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local3827dinitial_mpg40[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dgrav_acc_local3827dinitial_der121[((12 * 0) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 12, ddpotential_dx5670dinitial_inv_mpg71, ddpotential_dx5670dinitial_der181);
        }
      }
    }
    grav_acc_local122[0] = dpotential_dx182;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 1));
      if ((mappings_full_idx_symbol >= 36)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local3827dinitial_mpg40[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dgrav_acc_local3827dinitial_der121[((12 * 1) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 12, ddpotential_dy5697dinitial_inv_mpg27, ddpotential_dy5697dinitial_der184);
        }
      }
    }
    grav_acc_local122[1] = dpotential_dy185;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (12 * 2));
      if ((mappings_full_idx_symbol >= 36)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local3827dinitial_mpg40[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
          break;
        } else {
          dgrav_acc_local3827dinitial_der121[((12 * 2) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 12, ddpotential_dz5724dinitial_inv_mpg57, ddpotential_dz5724dinitial_der187);
        }
      }
    }
    grav_acc_local122[2] = dpotential_dz188;
    for (int k194 = 0; k194 < 3; k194++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * k194));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc3852dinitial_mpg68[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dgrav_acc3852dinitial_der123[((12 * k194) + mapped_idx)] = ((((get_dfdx_cell(0, 12, al_index_name_symbol, 12, dgrav_acc_local3827dinitial_inv_mpg41, dgrav_acc_local3827dinitial_der121) * rot127[((3 * k194) + 0)]) + (get_dfdx_cell(1, 12, al_index_name_symbol, 12, dgrav_acc_local3827dinitial_inv_mpg41, dgrav_acc_local3827dinitial_der121) * rot127[((3 * k194) + 1)])) + (get_dfdx_cell(2, 12, al_index_name_symbol, 12, dgrav_acc_local3827dinitial_inv_mpg41, dgrav_acc_local3827dinitial_der121) * rot127[((3 * k194) + 2)])) / 2.8432269415759565119e-08);
          }
        }
      }
      grav_acc124[k194] = ((((rot127[((3 * k194) + 0)] * grav_acc_local122[0]) + (rot127[((3 * k194) + 1)] * grav_acc_local122[1])) + (rot127[((3 * k194) + 2)] * grav_acc_local122[2])) / 2.8432269415759565119e-08);
    }
    for (int slice_idx = 0; slice_idx < (((i144 * 3) + 3) - (i144 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * ((i144 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2105dinitial_mpg78[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dsat_acc2105dinitial_der88[((12 * ((i144 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i144 * 3) + slice_idx), 12, al_index_name_symbol, 12, dsat_acc2105dinitial_inv_mpg79, dsat_acc2105dinitial_der88) + (get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dgrav_acc3852dinitial_inv_mpg69, dgrav_acc3852dinitial_der123) * 1.5240390322546755564e-08));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i144 * 3) + 3) - (i144 * 3)); slice_idx++) {
      sat_acc89[((i144 * 3) + slice_idx)] = (sat_acc89[((i144 * 3) + slice_idx)] + (1.5240390322546755564e-08 * grav_acc124[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2123dinitial_mpg58[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dcentral_acc2123dinitial_der90[((12 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dcentral_acc2123dinitial_inv_mpg59, dcentral_acc2123dinitial_der90) - (get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dgrav_acc3852dinitial_inv_mpg69, dgrav_acc3852dinitial_der123) * sat_gms25[i144]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc91[(0 + slice_idx)] = (central_acc91[(0 + slice_idx)] - (sat_gms25[i144] * grav_acc124[(0 + slice_idx)]));
    }
  }
  for (int i195 = 0; i195 < 2; i195++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (12 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 36)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc3871dinitial_mpg38[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 12))) {
            break;
          } else {
            dacc3871dinitial_der125[((12 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i195 * 3) + slice_idx), 12, al_index_name_symbol, 12, dsat_acc2105dinitial_inv_mpg79, dsat_acc2105dinitial_der88) - get_dfdx_cell((0 + slice_idx), 12, al_index_name_symbol, 12, dcentral_acc2123dinitial_inv_mpg59, dcentral_acc2123dinitial_der90));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc126[(0 + slice_idx)] = (sat_acc89[((i195 * 3) + slice_idx)] - central_acc91[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i195 * 6) + 3) - (i195 * 6)); slice_idx++) {
        nepsatsystem19[((i195 * 6) + slice_idx)] = state94[(((i195 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i195 * 6) + 6) - ((i195 * 6) + 3)); slice_idx++) {
        nepsatsystem19[(((i195 * 6) + 3) + slice_idx)] = acc126[(0 + slice_idx)];
      }
    }
    for (int j196 = 0; j196 < 3; j196++) {
      {
        for (int slice_idx = 0; slice_idx < ((12 + (((((i195 * 6) + j196) + 1) * 2) * 6)) - (12 + ((((i195 * 6) + j196) * 2) * 6))); slice_idx++) {
          nepsatsystem19[((12 + ((((i195 * 6) + j196) * 2) * 6)) + slice_idx)] = state_and_derivatives21[((12 + (((((i195 * 6) + j196) + 3) * 2) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((12 + (((((i195 * 6) + j196) + 4) * 2) * 6)) - (12 + (((((i195 * 6) + j196) + 3) * 2) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (12 + (((((i195 * 6) + j196) + 3) * 2) * 6)));
        nepsatsystem19[func_slice_idx] = get_dfdx_cell(j196, 12, al_index_name_symbol, 12, dacc3871dinitial_inv_mpg39, dacc3871dinitial_der125);
      }
    }
  }
  return 0;
}

int nepsatsystem_noderiv(long double *restrict nepsatsystem_noderiv197, long double t198, long double *restrict state199, long double *restrict central_pos200, long double *restrict perturb_gms201, long double *restrict perturb_pos202, long double *restrict sat_gms203) {
  long double sat_acc204[6] = { 0.0 };
  long double central_acc205[3] = { 0.0 };
  for (int i206 = 0; i206 < 2; i206++) {
    long double r207[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r207[(0 + slice_idx)] = state199[((i206 * 6) + slice_idx)];
      }
    }
    long double dist2208;
    dist2208 = (((r207[0] * r207[0]) + (r207[1] * r207[1])) + (r207[2] * r207[2]));
    long double dist3209;
    dist3209 = (dist2208 * sqrtl(dist2208));
    {
      for (int slice_idx = 0; slice_idx < (((i206 * 3) + 3) - (i206 * 3)); slice_idx++) {
        sat_acc204[((i206 * 3) + slice_idx)] = ((-1.5240390322546755564e-08 * r207[(0 + slice_idx)]) / dist3209);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc205[(0 + slice_idx)] = (central_acc205[(0 + slice_idx)] + ((sat_gms203[i206] * r207[(0 + slice_idx)]) / dist3209));
      }
    }
  }
  for (int i210 = 0; i210 < 4; i210++) {
    long double r211[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r211[(0 + slice_idx)] = (perturb_pos202[((i210 * 3) + slice_idx)] - central_pos200[(0 + slice_idx)]);
      }
    }
    long double dist2212;
    dist2212 = (((r211[0] * r211[0]) + (r211[1] * r211[1])) + (r211[2] * r211[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc205[(0 + slice_idx)] = (central_acc205[(0 + slice_idx)] + (((perturb_gms201[i210] * r211[(0 + slice_idx)]) / dist2212) / sqrtl(dist2212)));
      }
    }
    for (int j213 = 0; j213 < 2; j213++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r211[(0 + slice_idx)] = (perturb_pos202[((i210 * 3) + slice_idx)] - (state199[((j213 * 6) + slice_idx)] + central_pos200[(0 + slice_idx)]));
        }
      }
      dist2212 = (((r211[0] * r211[0]) + (r211[1] * r211[1])) + (r211[2] * r211[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j213 * 3) + 3) - (j213 * 3)); slice_idx++) {
          sat_acc204[((j213 * 3) + slice_idx)] = (sat_acc204[((j213 * 3) + slice_idx)] + (((perturb_gms201[i210] * r211[(0 + slice_idx)]) / dist2212) / sqrtl(dist2212)));
        }
      }
    }
  }
  for (int i214 = 1; i214 < 2; i214++) {
    long double r215[3] = { 0.0 };
    for (int j216 = 0; j216 < i214; j216++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r215[(0 + slice_idx)] = (state199[((j216 * 6) + slice_idx)] - state199[((i214 * 6) + slice_idx)]);
        }
      }
      long double dist2217;
      dist2217 = (((r215[0] * r215[0]) + (r215[1] * r215[1])) + (r215[2] * r215[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i214 * 3) + 3) - (i214 * 3)); slice_idx++) {
          sat_acc204[((i214 * 3) + slice_idx)] = (sat_acc204[((i214 * 3) + slice_idx)] + (((sat_gms203[j216] * r215[(0 + slice_idx)]) / dist2217) / sqrtl(dist2217)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j216 * 3) + 3) - (j216 * 3)); slice_idx++) {
          sat_acc204[((j216 * 3) + slice_idx)] = (sat_acc204[((j216 * 3) + slice_idx)] - (((sat_gms203[i214] * r215[(0 + slice_idx)]) / dist2217) / sqrtl(dist2217)));
        }
      }
    }
  }
  long double grav_acc_local218[3] = { 0.0 };
  long double grav_acc219[3] = { 0.0 };
  long double acc220[3] = { 0.0 };
  long double rot221[9] = { 0.0 };
  long double v_0222[5] = { 0.0 };
  long double dv_0_dx223[5] = { 0.0 };
  long double dv_0_dy224[5] = { 0.0 };
  long double dv_0_dz225[5] = { 0.0 };
  {
    long double neptune_rotation_matrix9095227[9] = { 0.0 };
    long double _inl_neptune_rotation_matrix_T229;
    _inl_neptune_rotation_matrix_T229 = (t198 / 36525.0e0);
    long double _inl_neptune_rotation_matrix_N230;
    _inl_neptune_rotation_matrix_N230 = (((357.85000000000002274e0 + (52.3160000000000025e0 * _inl_neptune_rotation_matrix_T229)) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_alpha_0231;
    _inl_neptune_rotation_matrix_alpha_0231 = (((299.36000000000001364e0 + (0.6999999999999999556e0 * sinl(_inl_neptune_rotation_matrix_N230))) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_delta_0232;
    _inl_neptune_rotation_matrix_delta_0232 = (((43.460000000000000853e0 - (0.5100000000000000089e0 * cosl(_inl_neptune_rotation_matrix_N230))) * 3.141592653589793116e0) / 180.0e0);
    long double _inl_neptune_rotation_matrix_W233;
    _inl_neptune_rotation_matrix_W233 = ((((249.97800000000000864e0 + (541.13977569999997286e0 * t198)) - (0.47999999999999998224e0 * sinl(_inl_neptune_rotation_matrix_N230))) * 3.141592653589793116e0) / 180.0e0);
    neptune_rotation_matrix9095227[0] = ((-sinl(_inl_neptune_rotation_matrix_alpha_0231) * cosl(_inl_neptune_rotation_matrix_W233)) - ((cosl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_delta_0232)) * sinl(_inl_neptune_rotation_matrix_W233)));
    neptune_rotation_matrix9095227[1] = ((sinl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_W233)) - ((cosl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_delta_0232)) * cosl(_inl_neptune_rotation_matrix_W233)));
    neptune_rotation_matrix9095227[2] = (cosl(_inl_neptune_rotation_matrix_alpha_0231) * cosl(_inl_neptune_rotation_matrix_delta_0232));
    neptune_rotation_matrix9095227[3] = ((cosl(_inl_neptune_rotation_matrix_alpha_0231) * cosl(_inl_neptune_rotation_matrix_W233)) - ((sinl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_delta_0232)) * sinl(_inl_neptune_rotation_matrix_W233)));
    neptune_rotation_matrix9095227[4] = ((-cosl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_W233)) - ((sinl(_inl_neptune_rotation_matrix_alpha_0231) * sinl(_inl_neptune_rotation_matrix_delta_0232)) * cosl(_inl_neptune_rotation_matrix_W233)));
    neptune_rotation_matrix9095227[5] = (cosl(_inl_neptune_rotation_matrix_delta_0232) * sinl(_inl_neptune_rotation_matrix_alpha_0231));
    neptune_rotation_matrix9095227[6] = (cosl(_inl_neptune_rotation_matrix_delta_0232) * sinl(_inl_neptune_rotation_matrix_W233));
    neptune_rotation_matrix9095227[7] = (cosl(_inl_neptune_rotation_matrix_delta_0232) * cosl(_inl_neptune_rotation_matrix_W233));
    neptune_rotation_matrix9095227[8] = sinl(_inl_neptune_rotation_matrix_delta_0232);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot221[(0 + slice_idx)] = neptune_rotation_matrix9095227[(0 + slice_idx)];
    }
  }
  for (int i234 = 0; i234 < 2; i234++) {
    long double x235;
    x235 = (state199[((i234 * 6) + 0)] / 0.00016861871015922155108e0);
    long double y236;
    y236 = (state199[((i234 * 6) + 1)] / 0.00016861871015922155108e0);
    long double z237;
    z237 = (state199[((i234 * 6) + 2)] / 0.00016861871015922155108e0);
    long double _x238;
    _x238 = (((rot221[0] * x235) + (rot221[3] * y236)) + (rot221[6] * z237));
    long double _y239;
    _y239 = (((rot221[1] * x235) + (rot221[4] * y236)) + (rot221[7] * z237));
    long double _z240;
    _z240 = (((rot221[2] * x235) + (rot221[5] * y236)) + (rot221[8] * z237));
    long double _r2241;
    _r2241 = (((_x238 * _x238) + (_y239 * _y239)) + (_z240 * _z240));
    long double _r242;
    _r242 = sqrtl(_r2241);
    long double _r3243;
    _r3243 = (_r242 * _r2241);
    long double _r4244;
    _r4244 = (_r2241 * _r2241);
    long double _r5245;
    _r5245 = (_r4244 * _r242);
    v_0222[0] = (1.0e0 / _r242);
    v_0222[1] = ((1.7320508075688772936e0 * _z240) / _r3243);
    dv_0_dx223[0] = (-_x238 / _r3243);
    dv_0_dy224[0] = (-_y239 / _r3243);
    dv_0_dz225[0] = (-_z240 / _r3243);
    dv_0_dx223[1] = (((-5.1961524227066318805e0 * _z240) * _x238) / _r5245);
    dv_0_dy224[1] = (((-5.1961524227066318805e0 * _z240) * _y239) / _r5245);
    dv_0_dz225[1] = (1.7320508075688772936e0 * ((1.0e0 / _r3243) - (((3.0e0 * _z240) * _z240) / _r5245)));
    for (int n246 = 2; n246 < 5; n246++) {
      long double coef1247;
      coef1247 = sqrtl(((((2.0e0 * ((long double) n246)) - 1.0e0) * ((long double) ((2 * n246) + 1))) / ((long double) (n246 * n246))));
      long double coef2248;
      coef2248 = sqrtl(((((((long double) n246) - 1.0e0) * ((long double) (n246 - 1))) * ((long double) ((2 * n246) + 1))) / ((long double) ((n246 * n246) * ((2 * n246) - 3)))));
      v_0222[n246] = ((((coef1247 * v_0222[(n246 - 1)]) * _z240) / _r2241) - ((coef2248 * v_0222[(n246 - 2)]) / _r2241));
      dv_0_dx223[n246] = (((coef1247 * _z240) * ((dv_0_dx223[(n246 - 1)] / _r2241) - (((v_0222[(n246 - 1)] * 2.0e0) * _x238) / _r4244))) - (coef2248 * ((dv_0_dx223[(n246 - 2)] / _r2241) - (((v_0222[(n246 - 2)] * 2.0e0) * _x238) / _r4244))));
      dv_0_dy224[n246] = (((coef1247 * _z240) * ((dv_0_dy224[(n246 - 1)] / _r2241) - (((v_0222[(n246 - 1)] * 2.0e0) * _y239) / _r4244))) - (coef2248 * ((dv_0_dy224[(n246 - 2)] / _r2241) - (((v_0222[(n246 - 2)] * 2.0e0) * _y239) / _r4244))));
      dv_0_dz225[n246] = ((coef1247 * (((dv_0_dz225[(n246 - 1)] * _z240) / _r2241) + (v_0222[(n246 - 1)] * ((1.0e0 / _r2241) - (((2.0e0 * _z240) * _z240) / _r4244))))) - (coef2248 * ((dv_0_dz225[(n246 - 2)] / _r2241) - (((v_0222[(n246 - 2)] * 2.0e0) * _z240) / _r4244))));
    }
    long double dpotential_dx249;
    dpotential_dx249 = ((long double) 0);
    long double dpotential_dy250;
    dpotential_dy250 = ((long double) 0);
    long double dpotential_dz251;
    dpotential_dz251 = ((long double) 0);
    for (int n252 = 2; n252 < 5; n252++) {
      dpotential_dx249 = (dpotential_dx249 + (dv_0_dx223[n252] * central_grav11[n252]));
      dpotential_dy250 = (dpotential_dy250 + (dv_0_dy224[n252] * central_grav11[n252]));
      dpotential_dz251 = (dpotential_dz251 + (dv_0_dz225[n252] * central_grav11[n252]));
    }
    grav_acc_local218[0] = dpotential_dx249;
    grav_acc_local218[1] = dpotential_dy250;
    grav_acc_local218[2] = dpotential_dz251;
    for (int k253 = 0; k253 < 3; k253++) {
      grav_acc219[k253] = ((((rot221[((3 * k253) + 0)] * grav_acc_local218[0]) + (rot221[((3 * k253) + 1)] * grav_acc_local218[1])) + (rot221[((3 * k253) + 2)] * grav_acc_local218[2])) / 2.8432269415759565119e-08);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i234 * 3) + 3) - (i234 * 3)); slice_idx++) {
        sat_acc204[((i234 * 3) + slice_idx)] = (sat_acc204[((i234 * 3) + slice_idx)] + (1.5240390322546755564e-08 * grav_acc219[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc205[(0 + slice_idx)] = (central_acc205[(0 + slice_idx)] - (sat_gms203[i234] * grav_acc219[(0 + slice_idx)]));
      }
    }
  }
  for (int i254 = 0; i254 < 2; i254++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc220[(0 + slice_idx)] = (sat_acc204[((i254 * 3) + slice_idx)] - central_acc205[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i254 * 6) + 3) - (i254 * 6)); slice_idx++) {
        nepsatsystem_noderiv197[((i254 * 6) + slice_idx)] = state199[(((i254 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i254 * 6) + 6) - ((i254 * 6) + 3)); slice_idx++) {
        nepsatsystem_noderiv197[(((i254 * 6) + 3) + slice_idx)] = acc220[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

