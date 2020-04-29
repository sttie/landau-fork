#include <math.h>

static inline long double get_dfdx_cell(int full_idx81545, int dx_mapped_size81546, int al_idx81547, int inv_mapping_period81548, int *restrict dx_idx_mappings81549, long double *restrict der_vec81550) {
  return ((al_idx81547 < inv_mapping_period81548) ? ((dx_idx_mappings81549[((inv_mapping_period81548 * full_idx81545) + al_idx81547)] >= 0) ? der_vec81550[((dx_mapped_size81546 * full_idx81545) + dx_idx_mappings81549[((inv_mapping_period81548 * full_idx81545) + al_idx81547)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx81545, int dx_mapped_size81546, int al_idx81547, int inv_mapping_period81548, int *restrict dx_idx_mappings81549, long double *restrict der_vec81550) {
  return ((al_idx81547 < inv_mapping_period81548) ? ((al_idx81547 == full_idx81545) ? 1.0 : ((dx_idx_mappings81549[((inv_mapping_period81548 * full_idx81545) + al_idx81547)] >= 0) ? der_vec81550[((dx_mapped_size81546 * full_idx81545) + dx_idx_mappings81549[((inv_mapping_period81548 * full_idx81545) + al_idx81547)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx81552, int inv_mapping_period81553, int *restrict dx_idx_mappings81554, long double *restrict der_vec81555) {
  return ((((al_idx81552 < inv_mapping_period81553) ? dx_idx_mappings81554[al_idx81552] : -1) >= 0) ? der_vec81555[((al_idx81552 < inv_mapping_period81553) ? dx_idx_mappings81554[al_idx81552] : -1)] : 0.0);
}

const long double pole_ra81556[2] = { 268.05659500000001572e0, -0.006498999999999999569e0 };
const long double pole_dec81557[2] = { 64.49530300000000693e0, 0.0024130000000000002142e0 };
const long double pm81558[2] = { 284.94999999999998863e0, 870.5359999999999445e0 };
const long double nut_prec_ra81559[5] = { 0.00011699999999999999788e0, 0.00093800000000000003167e0, 0.0014319999999999998962e0, 3.000000000000000076e-05, 0.0021500000000000000014e0 };
const long double nut_prec_dec81560[5] = { 5.0000000000000002396e-05, 0.00040400000000000000798e0, 0.0006170000000000000354e0, -1.29999999999999992e-05, 0.0009259999999999999568e0 };
const long double Jabcde_081561[5] = { 99.36071400000000153e0, 175.89536899999998809e0, 300.32316200000002482e0, 114.01230499999999779e0, 49.511251000000001454e0 };
const long double Jabcde_T81562[5] = { 4850.4045999999998457e0, 1191.9604999999999109e0, 262.54750000000001364e0, 6070.247599999999693e0, 64.29999999999999716e0 };
const long double central_grav81563[7] = { 0.0e0, 0.0e0, -0.0065724808672554692115e0, 0.0e0, 0.00019554099999999999462e0, 0.0e0, -9.4975767597683728686e-06 };

int jupiter_rotation_matrix(long double *restrict jupiter_rotation_matrix81564, long double t81565) {
  long double T81567;
  T81567 = (t81565 / 36525.0e0);
  long double alpha_081571;
  alpha_081571 = (268.05659500000001572e0 + (-0.006498999999999999569e0 * T81567));
  long double delta_081575;
  delta_081575 = (64.49530300000000693e0 + (0.0024130000000000002142e0 * T81567));
  long double W81579;
  W81579 = (((284.94999999999998863e0 + (870.5359999999999445e0 * t81565)) * 3.141592653589793116e0) / 180.0e0);
  for (int i81583 = 0; i81583 < 5; i81583++) {
    long double J81585;
    J81585 = (((Jabcde_081561[i81583] + (Jabcde_T81562[i81583] * T81567)) * 3.141592653589793116e0) / 180.0e0);
    alpha_081571 = (alpha_081571 + (nut_prec_ra81559[i81583] * sinl(J81585)));
    delta_081575 = (delta_081575 + (nut_prec_dec81560[i81583] * cosl(J81585)));
  }
  alpha_081571 = (alpha_081571 * 0.017453292519943295088e0);
  delta_081575 = (delta_081575 * 0.017453292519943295088e0);
  jupiter_rotation_matrix81564[0] = ((-sinl(alpha_081571) * cosl(W81579)) - ((cosl(alpha_081571) * sinl(delta_081575)) * sinl(W81579)));
  jupiter_rotation_matrix81564[1] = ((sinl(alpha_081571) * sinl(W81579)) - ((cosl(alpha_081571) * sinl(delta_081575)) * cosl(W81579)));
  jupiter_rotation_matrix81564[2] = (cosl(alpha_081571) * cosl(delta_081575));
  jupiter_rotation_matrix81564[3] = ((cosl(alpha_081571) * cosl(W81579)) - ((sinl(alpha_081571) * sinl(delta_081575)) * sinl(W81579)));
  jupiter_rotation_matrix81564[4] = ((-cosl(alpha_081571) * sinl(W81579)) - ((sinl(alpha_081571) * sinl(delta_081575)) * cosl(W81579)));
  jupiter_rotation_matrix81564[5] = (cosl(delta_081575) * sinl(alpha_081571));
  jupiter_rotation_matrix81564[6] = (cosl(delta_081575) * sinl(W81579));
  jupiter_rotation_matrix81564[7] = (cosl(delta_081575) * cosl(W81579));
  jupiter_rotation_matrix81564[8] = sinl(delta_081575);
  return 0;
}

int jupsatsystem(long double *restrict jupsatsystem81628, long double t81629, long double *restrict state_and_derivatives81630, long double *restrict central_pos81631, long double *restrict perturb_gms81632, long double *restrict perturb_pos81633, long double *restrict sat_gms81634) {
  static int d_r55216dinitial_mpg81636[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_inv_mpg81637[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_mpg81638[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_inv_mpg81639[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dstate2924dinitial_mpg81640[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2924dinitial_inv_mpg81641[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dx4794dinitial_mpg81642[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_inv_mpg81643[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_mpg81644[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_inv_mpg81645[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_mpg81646[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_inv_mpg81647[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_mpg81648[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_inv_mpg81649[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_mpg81650[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_inv_mpg81651[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_mpg81652[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_inv_mpg81653[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_mpg81654[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_inv_mpg81655[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_mpg81656[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_inv_mpg81657[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_mpg81658[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_inv_mpg81659[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_mpg81660[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_inv_mpg81661[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_mpg81662[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_inv_mpg81663[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_mpg81664[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_inv_mpg81665[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_mpg81666[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_inv_mpg81667[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_mpg81668[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_inv_mpg81669[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_mpg81670[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_inv_mpg81671[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_mpg81672[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_inv_mpg81673[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_mpg81674[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_inv_mpg81675[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_mpg81676[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_inv_mpg81677[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_mpg81678[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_inv_mpg81679[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_mpg81680[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_inv_mpg81681[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_mpg81682[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_inv_mpg81683[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_mpg81684[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_inv_mpg81685[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dv_04620dinitial_mpg81686[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_04620dinitial_inv_mpg81687[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int d_r45191dinitial_mpg81688[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_inv_mpg81689[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_mpg81690[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_inv_mpg81691[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_mpg81692[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_inv_mpg81693[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  long double sat_acc81695[12] = { 0.0 };
  long double dsat_acc2828dinitial_der81694[288] = { 0.0 };
  long double central_acc81697[3] = { 0.0 };
  long double dcentral_acc2846dinitial_der81696[72] = { 0.0 };
  long double state_derivatives_initial81698[576] = { 0.0 };
  long double state81700[24] = { 0.0 };
  long double dstate2924dinitial_der81699[576] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 576; slice_idx++) {
      state_derivatives_initial81698[(0 + slice_idx)] = state_and_derivatives81630[(24 + slice_idx)];
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
      if ((mappings_full_idx_symbol >= 576)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dstate2924dinitial_mpg81640[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dstate2924dinitial_der81699[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    state81700[(0 + slice_idx)] = state_and_derivatives81630[(0 + slice_idx)];
  }
  long double dist281708;
  long double ddist23094dinitial_der81707[24] = { 0.0 };
  long double dist381710;
  long double ddist33107dinitial_der81709[24] = { 0.0 };
  for (int i81711 = 0; i81711 < 24; i81711++) {
    for (int j81713 = 0; j81713 < 24; j81713++) {
      {
        int mapped_idx;
        mapped_idx = ((j81713 < 24) ? dstate2924dinitial_inv_mpg81641[((i81711 * 24) + j81713)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2924dinitial_der81699[((i81711 * 24) + mapped_idx)] = state_derivatives_initial81698[(((i81711 * 4) * 6) + j81713)];
        }
      }
    }
  }
  for (int i81716 = 0; i81716 < 4; i81716++) {
    long double r81719[3] = { 0.0 };
    long double dr3376dinitial_der81718[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3376dinitial_mpg81666[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3376dinitial_der81718[((24 * (0 + slice_idx)) + mapped_idx)] = get_dfdx_cell(((i81716 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81719[(0 + slice_idx)] = state81700[((i81716 * 6) + slice_idx)];
    }
    for (int loop_var81727 = 0; loop_var81727 < 24; loop_var81727++) {
      int mapped_idx;
      mapped_idx = loop_var81727;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81654[loop_var81727];
      ddist23094dinitial_der81707[mapped_idx] = ((((r81719[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718)) + (r81719[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718))) + ((r81719[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718)) + (r81719[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718)))) + ((r81719[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718)) + (r81719[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718))));
    }
    dist281708 = (((r81719[0] * r81719[0]) + (r81719[1] * r81719[1])) + (r81719[2] * r81719[2]));
    for (int loop_var81732 = 0; loop_var81732 < 24; loop_var81732++) {
      int mapped_idx;
      mapped_idx = loop_var81732;
      int al_index_name_symbol;
      al_index_name_symbol = ddist33107dinitial_mpg81646[loop_var81732];
      ddist33107dinitial_der81709[mapped_idx] = ((dist281708 * (0.5 * (sqrtl((1.0 / dist281708)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707)))) + (sqrtl(dist281708) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707)));
    }
    dist381710 = (dist281708 * sqrtl(dist281708));
    for (int slice_idx = 0; slice_idx < (((i81716 * 3) + 3) - (i81716 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81716 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81690[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81694[((24 * ((i81716 * 3) + slice_idx)) + mapped_idx)] = ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718) * -2.8247609439046209905e-07) * dist381710) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81647, ddist33107dinitial_der81709) * (r81719[(0 + slice_idx)] * -2.8247609439046209905e-07))) / (dist381710 * dist381710));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81716 * 3) + 3) - (i81716 * 3)); slice_idx++) {
      sat_acc81695[((i81716 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r81719[(0 + slice_idx)]) / dist381710);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81658[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81696[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81659, dcentral_acc2846dinitial_der81696) + ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81667, dr3376dinitial_der81718) * sat_gms81634[i81716]) * dist381710) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81647, ddist33107dinitial_der81709) * (r81719[(0 + slice_idx)] * sat_gms81634[i81716]))) / (dist381710 * dist381710)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81697[(0 + slice_idx)] = (central_acc81697[(0 + slice_idx)] + ((sat_gms81634[i81716] * r81719[(0 + slice_idx)]) / dist381710));
    }
  }
  for (int i81742 = 0; i81742 < 4; i81742++) {
    long double r81745[3] = { 0.0 };
    long double dr3708dinitial_der81744[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3708dinitial_mpg81668[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3708dinitial_der81744[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81745[(0 + slice_idx)] = (perturb_pos81633[((i81742 * 3) + slice_idx)] - central_pos81631[(0 + slice_idx)]);
    }
    for (int loop_var81752 = 0; loop_var81752 < 24; loop_var81752++) {
      int mapped_idx;
      mapped_idx = loop_var81752;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81654[loop_var81752];
      ddist23094dinitial_der81707[mapped_idx] = ((((r81745[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744))) + ((r81745[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)))) + ((r81745[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744))));
    }
    dist281708 = (((r81745[0] * r81745[0]) + (r81745[1] * r81745[1])) + (r81745[2] * r81745[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81658[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81696[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81659, dcentral_acc2846dinitial_der81696) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744) * perturb_gms81632[i81742]) * dist281708) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707) * (r81745[(0 + slice_idx)] * perturb_gms81632[i81742]))) / (dist281708 * dist281708)) * sqrtl(dist281708)) - ((0.5 * (sqrtl((1.0 / dist281708)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707))) * ((r81745[(0 + slice_idx)] * perturb_gms81632[i81742]) / dist281708))) / (sqrtl(dist281708) * sqrtl(dist281708))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81697[(0 + slice_idx)] = (central_acc81697[(0 + slice_idx)] + (((perturb_gms81632[i81742] * r81745[(0 + slice_idx)]) / dist281708) / sqrtl(dist281708)));
    }
    for (int j81758 = 0; j81758 < 4; j81758++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3708dinitial_mpg81668[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr3708dinitial_der81744[((24 * (0 + slice_idx)) + mapped_idx)] = -get_dfdx_cell(((j81758 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81745[(0 + slice_idx)] = (perturb_pos81633[((i81742 * 3) + slice_idx)] - (state81700[((j81758 * 6) + slice_idx)] + central_pos81631[(0 + slice_idx)]));
      }
      for (int loop_var81767 = 0; loop_var81767 < 24; loop_var81767++) {
        int mapped_idx;
        mapped_idx = loop_var81767;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81654[loop_var81767];
        ddist23094dinitial_der81707[mapped_idx] = ((((r81745[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744))) + ((r81745[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)))) + ((r81745[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744)) + (r81745[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744))));
      }
      dist281708 = (((r81745[0] * r81745[0]) + (r81745[1] * r81745[1])) + (r81745[2] * r81745[2]));
      for (int slice_idx = 0; slice_idx < (((j81758 * 3) + 3) - (j81758 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81758 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81690[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81694[((24 * ((j81758 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81758 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81691, dsat_acc2828dinitial_der81694) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81669, dr3708dinitial_der81744) * perturb_gms81632[i81742]) * dist281708) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707) * (r81745[(0 + slice_idx)] * perturb_gms81632[i81742]))) / (dist281708 * dist281708)) * sqrtl(dist281708)) - ((0.5 * (sqrtl((1.0 / dist281708)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707))) * ((r81745[(0 + slice_idx)] * perturb_gms81632[i81742]) / dist281708))) / (sqrtl(dist281708) * sqrtl(dist281708))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81758 * 3) + 3) - (j81758 * 3)); slice_idx++) {
        sat_acc81695[((j81758 * 3) + slice_idx)] = (sat_acc81695[((j81758 * 3) + slice_idx)] + (((perturb_gms81632[i81742] * r81745[(0 + slice_idx)]) / dist281708) / sqrtl(dist281708)));
      }
    }
  }
  for (int i81773 = 1; i81773 < 4; i81773++) {
    long double r81776[3] = { 0.0 };
    long double dr4219dinitial_der81775[72] = { 0.0 };
    for (int j81777 = 0; j81777 < i81773; j81777++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr4219dinitial_mpg81674[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr4219dinitial_der81775[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81777 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699) - get_dfdx_cell(((i81773 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81776[(0 + slice_idx)] = (state81700[((j81777 * 6) + slice_idx)] - state81700[((i81773 * 6) + slice_idx)]);
      }
      for (int loop_var81786 = 0; loop_var81786 < 24; loop_var81786++) {
        int mapped_idx;
        mapped_idx = loop_var81786;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81654[loop_var81786];
        ddist23094dinitial_der81707[mapped_idx] = ((((r81776[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775)) + (r81776[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775))) + ((r81776[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775)) + (r81776[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775)))) + ((r81776[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775)) + (r81776[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775))));
      }
      dist281708 = (((r81776[0] * r81776[0]) + (r81776[1] * r81776[1])) + (r81776[2] * r81776[2]));
      for (int slice_idx = 0; slice_idx < (((i81773 * 3) + 3) - (i81773 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((i81773 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81690[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81694[((24 * ((i81773 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81773 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81691, dsat_acc2828dinitial_der81694) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775) * sat_gms81634[j81777]) * dist281708) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707) * (r81776[(0 + slice_idx)] * sat_gms81634[j81777]))) / (dist281708 * dist281708)) * sqrtl(dist281708)) - ((0.5 * (sqrtl((1.0 / dist281708)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707))) * ((r81776[(0 + slice_idx)] * sat_gms81634[j81777]) / dist281708))) / (sqrtl(dist281708) * sqrtl(dist281708))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i81773 * 3) + 3) - (i81773 * 3)); slice_idx++) {
        sat_acc81695[((i81773 * 3) + slice_idx)] = (sat_acc81695[((i81773 * 3) + slice_idx)] + (((sat_gms81634[j81777] * r81776[(0 + slice_idx)]) / dist281708) / sqrtl(dist281708)));
      }
      for (int slice_idx = 0; slice_idx < (((j81777 * 3) + 3) - (j81777 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81777 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81690[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81694[((24 * ((j81777 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81777 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81691, dsat_acc2828dinitial_der81694) - (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81675, dr4219dinitial_der81775) * sat_gms81634[i81773]) * dist281708) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707) * (r81776[(0 + slice_idx)] * sat_gms81634[i81773]))) / (dist281708 * dist281708)) * sqrtl(dist281708)) - ((0.5 * (sqrtl((1.0 / dist281708)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81655, ddist23094dinitial_der81707))) * ((r81776[(0 + slice_idx)] * sat_gms81634[i81773]) / dist281708))) / (sqrtl(dist281708) * sqrtl(dist281708))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81777 * 3) + 3) - (j81777 * 3)); slice_idx++) {
        sat_acc81695[((j81777 * 3) + slice_idx)] = (sat_acc81695[((j81777 * 3) + slice_idx)] - (((sat_gms81634[i81773] * r81776[(0 + slice_idx)]) / dist281708) / sqrtl(dist281708)));
      }
    }
  }
  long double grav_acc_local81797[3] = { 0.0 };
  long double dgrav_acc_local4544dinitial_der81796[72] = { 0.0 };
  long double grav_acc81799[3] = { 0.0 };
  long double dgrav_acc4569dinitial_der81798[72] = { 0.0 };
  long double acc81801[3] = { 0.0 };
  long double dacc4588dinitial_der81800[72] = { 0.0 };
  long double rot81802[9] = { 0.0 };
  long double v_081804[7] = { 0.0 };
  long double dv_04620dinitial_der81803[168] = { 0.0 };
  long double dv_0_dx81806[7] = { 0.0 };
  long double ddv_0_dx4638dinitial_der81805[168] = { 0.0 };
  long double dv_0_dy81808[7] = { 0.0 };
  long double ddv_0_dy4660dinitial_der81807[168] = { 0.0 };
  long double dv_0_dz81810[7] = { 0.0 };
  long double ddv_0_dz4682dinitial_der81809[168] = { 0.0 };
  {
    long double jupiter_rotation_matrix81819[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix81819, t81629);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot81802[(0 + slice_idx)] = jupiter_rotation_matrix81819[slice_idx];
    }
  }
  for (int i81820 = 0; i81820 < 4; i81820++) {
    long double x81823;
    long double dx4794dinitial_der81822[24] = { 0.0 };
    for (int loop_var81827 = 0; loop_var81827 < 24; loop_var81827++) {
      int mapped_idx;
      mapped_idx = loop_var81827;
      int al_index_name_symbol;
      al_index_name_symbol = dx4794dinitial_mpg81642[loop_var81827];
      dx4794dinitial_der81822[mapped_idx] = (get_dfdx_cell(((i81820 * 6) + 0), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699) / 0.0004778945025452157572e0);
    }
    x81823 = (state81700[((i81820 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y81830;
    long double dy4843dinitial_der81829[24] = { 0.0 };
    for (int loop_var81834 = 0; loop_var81834 < 24; loop_var81834++) {
      int mapped_idx;
      mapped_idx = loop_var81834;
      int al_index_name_symbol;
      al_index_name_symbol = dy4843dinitial_mpg81670[loop_var81834];
      dy4843dinitial_der81829[mapped_idx] = (get_dfdx_cell(((i81820 * 6) + 1), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699) / 0.0004778945025452157572e0);
    }
    y81830 = (state81700[((i81820 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z81837;
    long double dz4892dinitial_der81836[24] = { 0.0 };
    for (int loop_var81841 = 0; loop_var81841 < 24; loop_var81841++) {
      int mapped_idx;
      mapped_idx = loop_var81841;
      int al_index_name_symbol;
      al_index_name_symbol = dz4892dinitial_mpg81664[loop_var81841];
      dz4892dinitial_der81836[mapped_idx] = (get_dfdx_cell(((i81820 * 6) + 2), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81641, dstate2924dinitial_der81699) / 0.0004778945025452157572e0);
    }
    z81837 = (state81700[((i81820 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x81844;
    long double d_x4942dinitial_der81843[24] = { 0.0 };
    for (int loop_var81848 = 0; loop_var81848 < 24; loop_var81848++) {
      int mapped_idx;
      mapped_idx = loop_var81848;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4942dinitial_mpg81676[loop_var81848];
      d_x4942dinitial_der81843[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81643, dx4794dinitial_der81822) * rot81802[0]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81671, dy4843dinitial_der81829) * rot81802[3])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81665, dz4892dinitial_der81836) * rot81802[6]));
    }
    _x81844 = (((rot81802[0] * x81823) + (rot81802[3] * y81830)) + (rot81802[6] * z81837));
    long double _y81851;
    long double d_y4993dinitial_der81850[24] = { 0.0 };
    for (int loop_var81855 = 0; loop_var81855 < 24; loop_var81855++) {
      int mapped_idx;
      mapped_idx = loop_var81855;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4993dinitial_mpg81644[loop_var81855];
      d_y4993dinitial_der81850[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81643, dx4794dinitial_der81822) * rot81802[1]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81671, dy4843dinitial_der81829) * rot81802[4])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81665, dz4892dinitial_der81836) * rot81802[7]));
    }
    _y81851 = (((rot81802[1] * x81823) + (rot81802[4] * y81830)) + (rot81802[7] * z81837));
    long double _z81858;
    long double d_z5044dinitial_der81857[24] = { 0.0 };
    for (int loop_var81862 = 0; loop_var81862 < 24; loop_var81862++) {
      int mapped_idx;
      mapped_idx = loop_var81862;
      int al_index_name_symbol;
      al_index_name_symbol = d_z5044dinitial_mpg81648[loop_var81862];
      d_z5044dinitial_der81857[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81643, dx4794dinitial_der81822) * rot81802[2]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81671, dy4843dinitial_der81829) * rot81802[5])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81665, dz4892dinitial_der81836) * rot81802[8]));
    }
    _z81858 = (((rot81802[2] * x81823) + (rot81802[5] * y81830)) + (rot81802[8] * z81837));
    long double _r281865;
    long double d_r25100dinitial_der81864[24] = { 0.0 };
    for (int loop_var81869 = 0; loop_var81869 < 24; loop_var81869++) {
      int mapped_idx;
      mapped_idx = loop_var81869;
      int al_index_name_symbol;
      al_index_name_symbol = d_r25100dinitial_mpg81660[loop_var81869];
      d_r25100dinitial_der81864[mapped_idx] = ((((_x81844 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843)) + (_x81844 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843))) + ((_y81851 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850)) + (_y81851 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850)))) + ((_z81858 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857))));
    }
    _r281865 = (((_x81844 * _x81844) + (_y81851 * _y81851)) + (_z81858 * _z81858));
    long double _r81872;
    long double d_r5143dinitial_der81871[24] = { 0.0 };
    for (int loop_var81876 = 0; loop_var81876 < 24; loop_var81876++) {
      int mapped_idx;
      mapped_idx = loop_var81876;
      int al_index_name_symbol;
      al_index_name_symbol = d_r5143dinitial_mpg81680[loop_var81876];
      d_r5143dinitial_der81871[mapped_idx] = (0.5 * (sqrtl((1.0 / _r281865)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864)));
    }
    _r81872 = sqrtl(_r281865);
    long double _r381879;
    long double d_r35167dinitial_der81878[24] = { 0.0 };
    for (int loop_var81883 = 0; loop_var81883 < 24; loop_var81883++) {
      int mapped_idx;
      mapped_idx = loop_var81883;
      int al_index_name_symbol;
      al_index_name_symbol = d_r35167dinitial_mpg81662[loop_var81883];
      d_r35167dinitial_der81878[mapped_idx] = ((_r81872 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864)) + (_r281865 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81681, d_r5143dinitial_der81871)));
    }
    _r381879 = (_r81872 * _r281865);
    long double _r481886;
    long double d_r45191dinitial_der81885[24] = { 0.0 };
    for (int loop_var81890 = 0; loop_var81890 < 24; loop_var81890++) {
      int mapped_idx;
      mapped_idx = loop_var81890;
      int al_index_name_symbol;
      al_index_name_symbol = d_r45191dinitial_mpg81688[loop_var81890];
      d_r45191dinitial_der81885[mapped_idx] = ((_r281865 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864)) + (_r281865 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864)));
    }
    _r481886 = (_r281865 * _r281865);
    long double _r581893;
    long double d_r55216dinitial_der81892[24] = { 0.0 };
    for (int loop_var81897 = 0; loop_var81897 < 24; loop_var81897++) {
      int mapped_idx;
      mapped_idx = loop_var81897;
      int al_index_name_symbol;
      al_index_name_symbol = d_r55216dinitial_mpg81636[loop_var81897];
      d_r55216dinitial_der81892[mapped_idx] = ((_r481886 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81681, d_r5143dinitial_der81871)) + (_r81872 * get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885)));
    }
    _r581893 = (_r481886 * _r81872);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81686[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81803[((24 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r81872 * _r81872)) * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81681, d_r5143dinitial_der81871))));
        }
      }
    }
    v_081804[0] = (1.0e0 / _r81872);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81686[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81803[((24 * 1) + mapped_idx)] = ((((get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * 1.7320508075688772936e0) * _r381879) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81663, d_r35167dinitial_der81878) * (_z81858 * 1.7320508075688772936e0))) / (_r381879 * _r381879));
        }
      }
    }
    v_081804[1] = ((1.7320508075688772936e0 * _z81858) / _r381879);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81650[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81805[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843) * _r381879) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81663, d_r35167dinitial_der81878) * -_x81844)) / (_r381879 * _r381879));
        }
      }
    }
    dv_0_dx81806[0] = (-_x81844 / _r381879);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81682[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81807[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850) * _r381879) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81663, d_r35167dinitial_der81878) * -_y81851)) / (_r381879 * _r381879));
        }
      }
    }
    dv_0_dy81808[0] = (-_y81851 / _r381879);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81656[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81809[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * _r381879) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81663, d_r35167dinitial_der81878) * -_z81858)) / (_r381879 * _r381879));
        }
      }
    }
    dv_0_dz81810[0] = (-_z81858 / _r381879);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81650[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81805[((24 * 1) + mapped_idx)] = ((((((_z81858 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843)) + (_x81844 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * -5.1961524227066318805e0))) * _r581893) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81637, d_r55216dinitial_der81892) * ((_z81858 * -5.1961524227066318805e0) * _x81844))) / (_r581893 * _r581893));
        }
      }
    }
    dv_0_dx81806[1] = (((-5.1961524227066318805e0 * _z81858) * _x81844) / _r581893);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81682[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81807[((24 * 1) + mapped_idx)] = ((((((_z81858 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850)) + (_y81851 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * -5.1961524227066318805e0))) * _r581893) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81637, d_r55216dinitial_der81892) * ((_z81858 * -5.1961524227066318805e0) * _y81851))) / (_r581893 * _r581893));
        }
      }
    }
    dv_0_dy81808[1] = (((-5.1961524227066318805e0 * _z81858) * _y81851) / _r581893);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81656[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81809[((24 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r381879 * _r381879)) * get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81663, d_r35167dinitial_der81878)))) - ((((((_z81858 * 3.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * 3.0e0))) * _r581893) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81637, d_r55216dinitial_der81892) * ((_z81858 * 3.0e0) * _z81858))) / (_r581893 * _r581893))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz81810[1] = (1.7320508075688772936e0 * ((1.0e0 / _r381879) - (((3.0e0 * _z81858) * _z81858) / _r581893)));
    for (int n81931 = 2; n81931 < 7; n81931++) {
      long double coef181933;
      coef181933 = sqrtl(((((2.0e0 * ((long double) n81931)) - 1.0e0) * ((long double) ((2 * n81931) + 1))) / ((long double) (n81931 * n81931))));
      long double coef281937;
      coef281937 = sqrtl(((((((long double) n81931) - 1.0e0) * ((long double) (n81931 - 1))) * ((long double) ((2 * n81931) + 1))) / ((long double) ((n81931 * n81931) * ((2 * n81931) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81931));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_04620dinitial_mpg81686[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dv_04620dinitial_der81803[((24 * n81931) + mapped_idx)] = (((((((v_081804[(n81931 - 1)] * coef181933) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * (get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * coef181933))) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * ((v_081804[(n81931 - 1)] * coef181933) * _z81858))) / (_r281865 * _r281865)) - ((((get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * coef281937) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * (v_081804[(n81931 - 2)] * coef281937))) / (_r281865 * _r281865)));
          }
        }
      }
      v_081804[n81931] = ((((coef181933 * v_081804[(n81931 - 1)]) * _z81858) / _r281865) - ((coef281937 * v_081804[(n81931 - 2)]) / _r281865));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81931));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx4638dinitial_mpg81650[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dx4638dinitial_der81805[((24 * n81931) + mapped_idx)] = ((((_z81858 * coef181933) * ((((get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81651, ddv_0_dx4638dinitial_der81805) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * dv_0_dx81806[(n81931 - 1)])) / (_r281865 * _r281865)) - ((((((v_081804[(n81931 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843)) + (_x81844 * (get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((v_081804[(n81931 - 1)] * 2.0e0) * _x81844))) / (_r481886 * _r481886)))) + (((dv_0_dx81806[(n81931 - 1)] / _r281865) - (((v_081804[(n81931 - 1)] * 2.0e0) * _x81844) / _r481886)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * coef181933))) - (((((get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81651, ddv_0_dx4638dinitial_der81805) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * dv_0_dx81806[(n81931 - 2)])) / (_r281865 * _r281865)) - ((((((v_081804[(n81931 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81677, d_x4942dinitial_der81843)) + (_x81844 * (get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((v_081804[(n81931 - 2)] * 2.0e0) * _x81844))) / (_r481886 * _r481886))) * coef281937));
          }
        }
      }
      dv_0_dx81806[n81931] = (((coef181933 * _z81858) * ((dv_0_dx81806[(n81931 - 1)] / _r281865) - (((v_081804[(n81931 - 1)] * 2.0e0) * _x81844) / _r481886))) - (coef281937 * ((dv_0_dx81806[(n81931 - 2)] / _r281865) - (((v_081804[(n81931 - 2)] * 2.0e0) * _x81844) / _r481886))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81931));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy4660dinitial_mpg81682[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dy4660dinitial_der81807[((24 * n81931) + mapped_idx)] = ((((_z81858 * coef181933) * ((((get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81683, ddv_0_dy4660dinitial_der81807) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * dv_0_dy81808[(n81931 - 1)])) / (_r281865 * _r281865)) - ((((((v_081804[(n81931 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850)) + (_y81851 * (get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((v_081804[(n81931 - 1)] * 2.0e0) * _y81851))) / (_r481886 * _r481886)))) + (((dv_0_dy81808[(n81931 - 1)] / _r281865) - (((v_081804[(n81931 - 1)] * 2.0e0) * _y81851) / _r481886)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * coef181933))) - (((((get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81683, ddv_0_dy4660dinitial_der81807) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * dv_0_dy81808[(n81931 - 2)])) / (_r281865 * _r281865)) - ((((((v_081804[(n81931 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81645, d_y4993dinitial_der81850)) + (_y81851 * (get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((v_081804[(n81931 - 2)] * 2.0e0) * _y81851))) / (_r481886 * _r481886))) * coef281937));
          }
        }
      }
      dv_0_dy81808[n81931] = (((coef181933 * _z81858) * ((dv_0_dy81808[(n81931 - 1)] / _r281865) - (((v_081804[(n81931 - 1)] * 2.0e0) * _y81851) / _r481886))) - (coef281937 * ((dv_0_dy81808[(n81931 - 2)] / _r281865) - (((v_081804[(n81931 - 2)] * 2.0e0) * _y81851) / _r481886))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81931));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz4682dinitial_mpg81656[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dz4682dinitial_der81809[((24 * n81931) + mapped_idx)] = ((((((((dv_0_dz81810[(n81931 - 1)] * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81657, ddv_0_dz4682dinitial_der81809))) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * (dv_0_dz81810[(n81931 - 1)] * _z81858))) / (_r281865 * _r281865)) + ((v_081804[(n81931 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r281865 * _r281865)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864)))) - ((((((_z81858 * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((_z81858 * 2.0e0) * _z81858))) / (_r481886 * _r481886)))) + (((1.0e0 / _r281865) - (((_z81858 * 2.0e0) * _z81858) / _r481886)) * get_dfdx_cell((n81931 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803)))) * coef181933) - (((((get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81657, ddv_0_dz4682dinitial_der81809) * _r281865) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81661, d_r25100dinitial_der81864) * dv_0_dz81810[(n81931 - 2)])) / (_r281865 * _r281865)) - ((((((v_081804[(n81931 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81649, d_z5044dinitial_der81857)) + (_z81858 * (get_dfdx_cell((n81931 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81687, dv_04620dinitial_der81803) * 2.0e0))) * _r481886) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81689, d_r45191dinitial_der81885) * ((v_081804[(n81931 - 2)] * 2.0e0) * _z81858))) / (_r481886 * _r481886))) * coef281937));
          }
        }
      }
      dv_0_dz81810[n81931] = ((coef181933 * (((dv_0_dz81810[(n81931 - 1)] * _z81858) / _r281865) + (v_081804[(n81931 - 1)] * ((1.0e0 / _r281865) - (((2.0e0 * _z81858) * _z81858) / _r481886))))) - (coef281937 * ((dv_0_dz81810[(n81931 - 2)] / _r281865) - (((v_081804[(n81931 - 2)] * 2.0e0) * _z81858) / _r481886))));
    }
    long double dpotential_dx81958;
    long double ddpotential_dx6387dinitial_der81957[24] = { 0.0 };
    for (int loop_var81962 = 0; loop_var81962 < 24; loop_var81962++) {
      int mapped_idx;
      mapped_idx = loop_var81962;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx6387dinitial_mpg81652[loop_var81962];
      ddpotential_dx6387dinitial_der81957[mapped_idx] = 0.0;
    }
    dpotential_dx81958 = ((long double) 0);
    long double dpotential_dy81964;
    long double ddpotential_dy6414dinitial_der81963[24] = { 0.0 };
    for (int loop_var81968 = 0; loop_var81968 < 24; loop_var81968++) {
      int mapped_idx;
      mapped_idx = loop_var81968;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy6414dinitial_mpg81692[loop_var81968];
      ddpotential_dy6414dinitial_der81963[mapped_idx] = 0.0;
    }
    dpotential_dy81964 = ((long double) 0);
    long double dpotential_dz81970;
    long double ddpotential_dz6441dinitial_der81969[24] = { 0.0 };
    for (int loop_var81974 = 0; loop_var81974 < 24; loop_var81974++) {
      int mapped_idx;
      mapped_idx = loop_var81974;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz6441dinitial_mpg81684[loop_var81974];
      ddpotential_dz6441dinitial_der81969[mapped_idx] = 0.0;
    }
    dpotential_dz81970 = ((long double) 0);
    for (int n81975 = 2; n81975 < 7; n81975++) {
      for (int loop_var81980 = 0; loop_var81980 < 24; loop_var81980++) {
        int mapped_idx;
        mapped_idx = loop_var81980;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx6387dinitial_mpg81652[loop_var81980];
        ddpotential_dx6387dinitial_der81957[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81653, ddpotential_dx6387dinitial_der81957) + (get_dfdx_cell(n81975, 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81651, ddv_0_dx4638dinitial_der81805) * central_grav81563[n81975]));
      }
      dpotential_dx81958 = (dpotential_dx81958 + (dv_0_dx81806[n81975] * central_grav81563[n81975]));
      for (int loop_var81985 = 0; loop_var81985 < 24; loop_var81985++) {
        int mapped_idx;
        mapped_idx = loop_var81985;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy6414dinitial_mpg81692[loop_var81985];
        ddpotential_dy6414dinitial_der81963[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81693, ddpotential_dy6414dinitial_der81963) + (get_dfdx_cell(n81975, 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81683, ddv_0_dy4660dinitial_der81807) * central_grav81563[n81975]));
      }
      dpotential_dy81964 = (dpotential_dy81964 + (dv_0_dy81808[n81975] * central_grav81563[n81975]));
      for (int loop_var81990 = 0; loop_var81990 < 24; loop_var81990++) {
        int mapped_idx;
        mapped_idx = loop_var81990;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz6441dinitial_mpg81684[loop_var81990];
        ddpotential_dz6441dinitial_der81969[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81685, ddpotential_dz6441dinitial_der81969) + (get_dfdx_cell(n81975, 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81657, ddv_0_dz4682dinitial_der81809) * central_grav81563[n81975]));
      }
      dpotential_dz81970 = (dpotential_dz81970 + (dv_0_dz81810[n81975] * central_grav81563[n81975]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81672[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81796[((24 * 0) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81653, ddpotential_dx6387dinitial_der81957);
        }
      }
    }
    grav_acc_local81797[0] = dpotential_dx81958;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81672[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81796[((24 * 1) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81693, ddpotential_dy6414dinitial_der81963);
        }
      }
    }
    grav_acc_local81797[1] = dpotential_dy81964;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 2));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81672[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81796[((24 * 2) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81685, ddpotential_dz6441dinitial_der81969);
        }
      }
    }
    grav_acc_local81797[2] = dpotential_dz81970;
    for (int k82004 = 0; k82004 < 3; k82004++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * k82004));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc4569dinitial_mpg81678[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dgrav_acc4569dinitial_der81798[((24 * k82004) + mapped_idx)] = ((((get_dfdx_cell(0, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81673, dgrav_acc_local4544dinitial_der81796) * rot81802[((3 * k82004) + 0)]) + (get_dfdx_cell(1, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81673, dgrav_acc_local4544dinitial_der81796) * rot81802[((3 * k82004) + 1)])) + (get_dfdx_cell(2, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81673, dgrav_acc_local4544dinitial_der81796) * rot81802[((3 * k82004) + 2)])) / 2.2838315556293922983e-07);
          }
        }
      }
      grav_acc81799[k82004] = ((((rot81802[((3 * k82004) + 0)] * grav_acc_local81797[0]) + (rot81802[((3 * k82004) + 1)] * grav_acc_local81797[1])) + (rot81802[((3 * k82004) + 2)] * grav_acc_local81797[2])) / 2.2838315556293922983e-07);
    }
    for (int slice_idx = 0; slice_idx < (((i81820 * 3) + 3) - (i81820 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81820 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81690[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81694[((24 * ((i81820 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81820 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81691, dsat_acc2828dinitial_der81694) + (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81679, dgrav_acc4569dinitial_der81798) * 2.8247609439046209905e-07));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81820 * 3) + 3) - (i81820 * 3)); slice_idx++) {
      sat_acc81695[((i81820 * 3) + slice_idx)] = (sat_acc81695[((i81820 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc81799[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81658[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81696[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81659, dcentral_acc2846dinitial_der81696) - (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81679, dgrav_acc4569dinitial_der81798) * sat_gms81634[i81820]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81697[(0 + slice_idx)] = (central_acc81697[(0 + slice_idx)] - (sat_gms81634[i81820] * grav_acc81799[(0 + slice_idx)]));
    }
  }
  for (int i82018 = 0; i82018 < 4; i82018++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc4588dinitial_mpg81638[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dacc4588dinitial_der81800[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82018 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81691, dsat_acc2828dinitial_der81694) - get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81659, dcentral_acc2846dinitial_der81696));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc81801[(0 + slice_idx)] = (sat_acc81695[((i82018 * 3) + slice_idx)] - central_acc81697[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82018 * 6) + 3) - (i82018 * 6)); slice_idx++) {
        jupsatsystem81628[((i82018 * 6) + slice_idx)] = state81700[(((i82018 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82018 * 6) + 6) - ((i82018 * 6) + 3)); slice_idx++) {
        jupsatsystem81628[(((i82018 * 6) + 3) + slice_idx)] = acc81801[(0 + slice_idx)];
      }
    }
    for (int j82030 = 0; j82030 < 3; j82030++) {
      {
        for (int slice_idx = 0; slice_idx < ((24 + (((((i82018 * 6) + j82030) + 1) * 4) * 6)) - (24 + ((((i82018 * 6) + j82030) * 4) * 6))); slice_idx++) {
          jupsatsystem81628[((24 + ((((i82018 * 6) + j82030) * 4) * 6)) + slice_idx)] = state_and_derivatives81630[((24 + (((((i82018 * 6) + j82030) + 3) * 4) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((24 + (((((i82018 * 6) + j82030) + 4) * 4) * 6)) - (24 + (((((i82018 * 6) + j82030) + 3) * 4) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (24 + (((((i82018 * 6) + j82030) + 3) * 4) * 6)));
        jupsatsystem81628[func_slice_idx] = get_dfdx_cell(j82030, 24, al_index_name_symbol, 24, dacc4588dinitial_inv_mpg81639, dacc4588dinitial_der81800);
      }
    }
  }
  return 0;
}

int jupsatsystem_noderiv(long double *restrict jupsatsystem_noderiv82037, long double t82038, long double *restrict state82039, long double *restrict central_pos82040, long double *restrict perturb_gms82041, long double *restrict perturb_pos82042, long double *restrict sat_gms82043) {
  long double sat_acc82045[12] = { 0.0 };
  long double central_acc82046[3] = { 0.0 };
  long double dist282047;
  long double dist382048;
  for (int i82049 = 0; i82049 < 4; i82049++) {
    long double r82051[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82051[(0 + slice_idx)] = state82039[((i82049 * 6) + slice_idx)];
      }
    }
    dist282047 = (((r82051[0] * r82051[0]) + (r82051[1] * r82051[1])) + (r82051[2] * r82051[2]));
    dist382048 = (dist282047 * sqrtl(dist282047));
    {
      for (int slice_idx = 0; slice_idx < (((i82049 * 3) + 3) - (i82049 * 3)); slice_idx++) {
        sat_acc82045[((i82049 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r82051[(0 + slice_idx)]) / dist382048);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82046[(0 + slice_idx)] = (central_acc82046[(0 + slice_idx)] + ((sat_gms82043[i82049] * r82051[(0 + slice_idx)]) / dist382048));
      }
    }
  }
  for (int i82067 = 0; i82067 < 4; i82067++) {
    long double r82069[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82069[(0 + slice_idx)] = (perturb_pos82042[((i82067 * 3) + slice_idx)] - central_pos82040[(0 + slice_idx)]);
      }
    }
    dist282047 = (((r82069[0] * r82069[0]) + (r82069[1] * r82069[1])) + (r82069[2] * r82069[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82046[(0 + slice_idx)] = (central_acc82046[(0 + slice_idx)] + (((perturb_gms82041[i82067] * r82069[(0 + slice_idx)]) / dist282047) / sqrtl(dist282047)));
      }
    }
    for (int j82079 = 0; j82079 < 4; j82079++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82069[(0 + slice_idx)] = (perturb_pos82042[((i82067 * 3) + slice_idx)] - (state82039[((j82079 * 6) + slice_idx)] + central_pos82040[(0 + slice_idx)]));
        }
      }
      dist282047 = (((r82069[0] * r82069[0]) + (r82069[1] * r82069[1])) + (r82069[2] * r82069[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j82079 * 3) + 3) - (j82079 * 3)); slice_idx++) {
          sat_acc82045[((j82079 * 3) + slice_idx)] = (sat_acc82045[((j82079 * 3) + slice_idx)] + (((perturb_gms82041[i82067] * r82069[(0 + slice_idx)]) / dist282047) / sqrtl(dist282047)));
        }
      }
    }
  }
  for (int i82090 = 1; i82090 < 4; i82090++) {
    long double r82092[3] = { 0.0 };
    for (int j82093 = 0; j82093 < i82090; j82093++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82092[(0 + slice_idx)] = (state82039[((j82093 * 6) + slice_idx)] - state82039[((i82090 * 6) + slice_idx)]);
        }
      }
      dist282047 = (((r82092[0] * r82092[0]) + (r82092[1] * r82092[1])) + (r82092[2] * r82092[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i82090 * 3) + 3) - (i82090 * 3)); slice_idx++) {
          sat_acc82045[((i82090 * 3) + slice_idx)] = (sat_acc82045[((i82090 * 3) + slice_idx)] + (((sat_gms82043[j82093] * r82092[(0 + slice_idx)]) / dist282047) / sqrtl(dist282047)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j82093 * 3) + 3) - (j82093 * 3)); slice_idx++) {
          sat_acc82045[((j82093 * 3) + slice_idx)] = (sat_acc82045[((j82093 * 3) + slice_idx)] - (((sat_gms82043[i82090] * r82092[(0 + slice_idx)]) / dist282047) / sqrtl(dist282047)));
        }
      }
    }
  }
  long double grav_acc_local82107[3] = { 0.0 };
  long double grav_acc82108[3] = { 0.0 };
  long double acc82109[3] = { 0.0 };
  long double rot82110[9] = { 0.0 };
  long double v_082111[7] = { 0.0 };
  long double dv_0_dx82112[7] = { 0.0 };
  long double dv_0_dy82113[7] = { 0.0 };
  long double dv_0_dz82114[7] = { 0.0 };
  {
    long double jupiter_rotation_matrix82123[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix82123, t82038);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot82110[(0 + slice_idx)] = jupiter_rotation_matrix82123[slice_idx];
    }
  }
  for (int i82124 = 0; i82124 < 4; i82124++) {
    long double x82126;
    x82126 = (state82039[((i82124 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y82130;
    y82130 = (state82039[((i82124 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z82134;
    z82134 = (state82039[((i82124 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82138;
    _x82138 = (((rot82110[0] * x82126) + (rot82110[3] * y82130)) + (rot82110[6] * z82134));
    long double _y82142;
    _y82142 = (((rot82110[1] * x82126) + (rot82110[4] * y82130)) + (rot82110[7] * z82134));
    long double _z82146;
    _z82146 = (((rot82110[2] * x82126) + (rot82110[5] * y82130)) + (rot82110[8] * z82134));
    long double _r282150;
    _r282150 = (((_x82138 * _x82138) + (_y82142 * _y82142)) + (_z82146 * _z82146));
    long double r82154;
    r82154 = sqrtl(_r282150);
    long double _r382158;
    _r382158 = (r82154 * _r282150);
    long double _r482162;
    _r482162 = (_r282150 * _r282150);
    long double _r582166;
    _r582166 = (_r482162 * r82154);
    v_082111[0] = (1.0e0 / r82154);
    v_082111[1] = ((1.7320508075688772936e0 * _z82146) / _r382158);
    dv_0_dx82112[0] = (-_x82138 / _r382158);
    dv_0_dy82113[0] = (-_y82142 / _r382158);
    dv_0_dz82114[0] = (-_z82146 / _r382158);
    dv_0_dx82112[1] = (((-5.1961524227066318805e0 * _z82146) * _x82138) / _r582166);
    dv_0_dy82113[1] = (((-5.1961524227066318805e0 * _z82146) * _y82142) / _r582166);
    dv_0_dz82114[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382158) - (((3.0e0 * _z82146) * _z82146) / _r582166)));
    for (int n82194 = 2; n82194 < 7; n82194++) {
      long double coef182196;
      coef182196 = sqrtl(((((2.0e0 * ((long double) n82194)) - 1.0e0) * ((long double) ((2 * n82194) + 1))) / ((long double) (n82194 * n82194))));
      long double coef282200;
      coef282200 = sqrtl(((((((long double) n82194) - 1.0e0) * ((long double) (n82194 - 1))) * ((long double) ((2 * n82194) + 1))) / ((long double) ((n82194 * n82194) * ((2 * n82194) - 3)))));
      v_082111[n82194] = ((((coef182196 * v_082111[(n82194 - 1)]) * _z82146) / _r282150) - ((coef282200 * v_082111[(n82194 - 2)]) / _r282150));
      dv_0_dx82112[n82194] = (((coef182196 * _z82146) * ((dv_0_dx82112[(n82194 - 1)] / _r282150) - (((v_082111[(n82194 - 1)] * 2.0e0) * _x82138) / _r482162))) - (coef282200 * ((dv_0_dx82112[(n82194 - 2)] / _r282150) - (((v_082111[(n82194 - 2)] * 2.0e0) * _x82138) / _r482162))));
      dv_0_dy82113[n82194] = (((coef182196 * _z82146) * ((dv_0_dy82113[(n82194 - 1)] / _r282150) - (((v_082111[(n82194 - 1)] * 2.0e0) * _y82142) / _r482162))) - (coef282200 * ((dv_0_dy82113[(n82194 - 2)] / _r282150) - (((v_082111[(n82194 - 2)] * 2.0e0) * _y82142) / _r482162))));
      dv_0_dz82114[n82194] = ((coef182196 * (((dv_0_dz82114[(n82194 - 1)] * _z82146) / _r282150) + (v_082111[(n82194 - 1)] * ((1.0e0 / _r282150) - (((2.0e0 * _z82146) * _z82146) / _r482162))))) - (coef282200 * ((dv_0_dz82114[(n82194 - 2)] / _r282150) - (((v_082111[(n82194 - 2)] * 2.0e0) * _z82146) / _r482162))));
    }
    long double dpotential_dx82216;
    dpotential_dx82216 = ((long double) 0);
    long double dpotential_dy82220;
    dpotential_dy82220 = ((long double) 0);
    long double dpotential_dz82224;
    dpotential_dz82224 = ((long double) 0);
    for (int n82228 = 2; n82228 < 7; n82228++) {
      dpotential_dx82216 = (dpotential_dx82216 + (dv_0_dx82112[n82228] * central_grav81563[n82228]));
      dpotential_dy82220 = (dpotential_dy82220 + (dv_0_dy82113[n82228] * central_grav81563[n82228]));
      dpotential_dz82224 = (dpotential_dz82224 + (dv_0_dz82114[n82228] * central_grav81563[n82228]));
    }
    grav_acc_local82107[0] = dpotential_dx82216;
    grav_acc_local82107[1] = dpotential_dy82220;
    grav_acc_local82107[2] = dpotential_dz82224;
    for (int k82248 = 0; k82248 < 3; k82248++) {
      grav_acc82108[k82248] = ((((rot82110[((3 * k82248) + 0)] * grav_acc_local82107[0]) + (rot82110[((3 * k82248) + 1)] * grav_acc_local82107[1])) + (rot82110[((3 * k82248) + 2)] * grav_acc_local82107[2])) / 2.2838315556293922983e-07);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82124 * 3) + 3) - (i82124 * 3)); slice_idx++) {
        sat_acc82045[((i82124 * 3) + slice_idx)] = (sat_acc82045[((i82124 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc82108[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82046[(0 + slice_idx)] = (central_acc82046[(0 + slice_idx)] - (sat_gms82043[i82124] * grav_acc82108[(0 + slice_idx)]));
      }
    }
  }
  for (int i82259 = 0; i82259 < 4; i82259++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc82109[(0 + slice_idx)] = (sat_acc82045[((i82259 * 3) + slice_idx)] - central_acc82046[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82259 * 6) + 3) - (i82259 * 6)); slice_idx++) {
        jupsatsystem_noderiv82037[((i82259 * 6) + slice_idx)] = state82039[(((i82259 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82259 * 6) + 6) - ((i82259 * 6) + 3)); slice_idx++) {
        jupsatsystem_noderiv82037[(((i82259 * 6) + 3) + slice_idx)] = acc82109[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

