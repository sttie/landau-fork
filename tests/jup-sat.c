#include <math.h>

static inline long double get_dfdx_cell(int full_idx81704, int dx_mapped_size81705, int al_idx81706, int inv_mapping_period81707, int *restrict dx_idx_mappings81708, long double *restrict der_vec81709) {
  return ((al_idx81706 < inv_mapping_period81707) ? ((dx_idx_mappings81708[((inv_mapping_period81707 * full_idx81704) + al_idx81706)] >= 0) ? der_vec81709[((dx_mapped_size81705 * full_idx81704) + dx_idx_mappings81708[((inv_mapping_period81707 * full_idx81704) + al_idx81706)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx81704, int dx_mapped_size81705, int al_idx81706, int inv_mapping_period81707, int *restrict dx_idx_mappings81708, long double *restrict der_vec81709) {
  return ((al_idx81706 < inv_mapping_period81707) ? ((al_idx81706 == full_idx81704) ? 1.0 : ((dx_idx_mappings81708[((inv_mapping_period81707 * full_idx81704) + al_idx81706)] >= 0) ? der_vec81709[((dx_mapped_size81705 * full_idx81704) + dx_idx_mappings81708[((inv_mapping_period81707 * full_idx81704) + al_idx81706)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx81711, int inv_mapping_period81712, int *restrict dx_idx_mappings81713, long double *restrict der_vec81714) {
  return ((((al_idx81711 < inv_mapping_period81712) ? dx_idx_mappings81713[al_idx81711] : -1) >= 0) ? der_vec81714[((al_idx81711 < inv_mapping_period81712) ? dx_idx_mappings81713[al_idx81711] : -1)] : 0.0);
}

const long double pole_ra81715[2] = { 268.05659500000001572e0, -0.006498999999999999569e0 };
const long double pole_dec81716[2] = { 64.49530300000000693e0, 0.0024130000000000002142e0 };
const long double pm81717[2] = { 284.94999999999998863e0, 870.5359999999999445e0 };
const long double nut_prec_ra81718[5] = { 0.00011699999999999999788e0, 0.00093800000000000003167e0, 0.0014319999999999998962e0, 3.000000000000000076e-05, 0.0021500000000000000014e0 };
const long double nut_prec_dec81719[5] = { 5.0000000000000002396e-05, 0.00040400000000000000798e0, 0.0006170000000000000354e0, -1.29999999999999992e-05, 0.0009259999999999999568e0 };
const long double Jabcde_081720[5] = { 99.36071400000000153e0, 175.89536899999998809e0, 300.32316200000002482e0, 114.01230499999999779e0, 49.511251000000001454e0 };
const long double Jabcde_T81721[5] = { 4850.4045999999998457e0, 1191.9604999999999109e0, 262.54750000000001364e0, 6070.247599999999693e0, 64.29999999999999716e0 };
const long double central_grav81722[7] = { 0.0e0, 0.0e0, -0.0065724808672554692115e0, 0.0e0, 0.00019554099999999999462e0, 0.0e0, -9.4975767597683728686e-06 };

int jupiter_rotation_matrix(long double *restrict jupiter_rotation_matrix81723, long double t81724) {
  long double T81726;
  T81726 = (t81724 / 36525.0e0);
  long double alpha_081730;
  alpha_081730 = (268.05659500000001572e0 + (-0.006498999999999999569e0 * T81726));
  long double delta_081734;
  delta_081734 = (64.49530300000000693e0 + (0.0024130000000000002142e0 * T81726));
  long double W81738;
  W81738 = (((284.94999999999998863e0 + (870.5359999999999445e0 * t81724)) * 3.141592653589793116e0) / 180.0e0);
  for (int i81742 = 0; i81742 < 5; i81742++) {
    long double J81744;
    J81744 = (((Jabcde_081720[i81742] + (Jabcde_T81721[i81742] * T81726)) * 3.141592653589793116e0) / 180.0e0);
    alpha_081730 = (alpha_081730 + (nut_prec_ra81718[i81742] * sinl(J81744)));
    delta_081734 = (delta_081734 + (nut_prec_dec81719[i81742] * cosl(J81744)));
  }
  alpha_081730 = (alpha_081730 * 0.017453292519943295088e0);
  delta_081734 = (delta_081734 * 0.017453292519943295088e0);
  jupiter_rotation_matrix81723[0] = ((-sinl(alpha_081730) * cosl(W81738)) - ((cosl(alpha_081730) * sinl(delta_081734)) * sinl(W81738)));
  jupiter_rotation_matrix81723[1] = ((sinl(alpha_081730) * sinl(W81738)) - ((cosl(alpha_081730) * sinl(delta_081734)) * cosl(W81738)));
  jupiter_rotation_matrix81723[2] = (cosl(alpha_081730) * cosl(delta_081734));
  jupiter_rotation_matrix81723[3] = ((cosl(alpha_081730) * cosl(W81738)) - ((sinl(alpha_081730) * sinl(delta_081734)) * sinl(W81738)));
  jupiter_rotation_matrix81723[4] = ((-cosl(alpha_081730) * sinl(W81738)) - ((sinl(alpha_081730) * sinl(delta_081734)) * cosl(W81738)));
  jupiter_rotation_matrix81723[5] = (cosl(delta_081734) * sinl(alpha_081730));
  jupiter_rotation_matrix81723[6] = (cosl(delta_081734) * sinl(W81738));
  jupiter_rotation_matrix81723[7] = (cosl(delta_081734) * cosl(W81738));
  jupiter_rotation_matrix81723[8] = sinl(delta_081734);
  return 0;
}

int jupsatsystem(long double *restrict jupsatsystem81787, long double t81788, long double *restrict state_and_derivatives81789, long double *restrict central_pos81790, long double *restrict perturb_gms81791, long double *restrict perturb_pos81792, long double *restrict sat_gms81793) {
  static int d_r55216dinitial_mpg81795[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_inv_mpg81796[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_mpg81797[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_inv_mpg81798[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dstate2924dinitial_mpg81799[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2924dinitial_inv_mpg81800[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dx4794dinitial_mpg81801[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_inv_mpg81802[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_mpg81803[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_inv_mpg81804[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_mpg81805[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_inv_mpg81806[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_mpg81807[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_inv_mpg81808[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_mpg81809[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_inv_mpg81810[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_mpg81811[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_inv_mpg81812[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_mpg81813[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_inv_mpg81814[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_mpg81815[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_inv_mpg81816[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_mpg81817[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_inv_mpg81818[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_mpg81819[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_inv_mpg81820[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_mpg81821[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_inv_mpg81822[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_mpg81823[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_inv_mpg81824[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_mpg81825[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_inv_mpg81826[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_mpg81827[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_inv_mpg81828[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_mpg81829[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_inv_mpg81830[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_mpg81831[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_inv_mpg81832[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_mpg81833[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_inv_mpg81834[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_mpg81835[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_inv_mpg81836[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_mpg81837[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_inv_mpg81838[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_mpg81839[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_inv_mpg81840[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_mpg81841[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_inv_mpg81842[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_mpg81843[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_inv_mpg81844[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dv_04620dinitial_mpg81845[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_04620dinitial_inv_mpg81846[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int d_r45191dinitial_mpg81847[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_inv_mpg81848[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_mpg81849[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_inv_mpg81850[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_mpg81851[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_inv_mpg81852[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  long double sat_acc81854[12] = { 0.0 };
  long double dsat_acc2828dinitial_der81853[288] = { 0.0 };
  long double central_acc81856[3] = { 0.0 };
  long double dcentral_acc2846dinitial_der81855[72] = { 0.0 };
  long double state_derivatives_initial81857[576] = { 0.0 };
  long double state81859[24] = { 0.0 };
  long double dstate2924dinitial_der81858[576] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 576; slice_idx++) {
      state_derivatives_initial81857[(0 + slice_idx)] = state_and_derivatives81789[(24 + slice_idx)];
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
        al_index_name_symbol = dstate2924dinitial_mpg81799[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dstate2924dinitial_der81858[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    state81859[(0 + slice_idx)] = state_and_derivatives81789[(0 + slice_idx)];
  }
  long double dist281867;
  long double ddist23094dinitial_der81866[24] = { 0.0 };
  long double dist381869;
  long double ddist33107dinitial_der81868[24] = { 0.0 };
  for (int i81870 = 0; i81870 < 24; i81870++) {
    for (int j81872 = 0; j81872 < 24; j81872++) {
      {
        int mapped_idx;
        mapped_idx = ((j81872 < 24) ? dstate2924dinitial_inv_mpg81800[((i81870 * 24) + j81872)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2924dinitial_der81858[((i81870 * 24) + mapped_idx)] = state_derivatives_initial81857[(((i81870 * 4) * 6) + j81872)];
        }
      }
    }
  }
  for (int i81875 = 0; i81875 < 4; i81875++) {
    long double r81878[3] = { 0.0 };
    long double dr3376dinitial_der81877[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3376dinitial_mpg81825[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3376dinitial_der81877[((24 * (0 + slice_idx)) + mapped_idx)] = get_dfdx_cell(((i81875 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81878[(0 + slice_idx)] = state81859[((i81875 * 6) + slice_idx)];
    }
    for (int loop_var81886 = 0; loop_var81886 < 24; loop_var81886++) {
      int mapped_idx;
      mapped_idx = loop_var81886;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81813[loop_var81886];
      ddist23094dinitial_der81866[mapped_idx] = ((((r81878[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877)) + (r81878[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877))) + ((r81878[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877)) + (r81878[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877)))) + ((r81878[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877)) + (r81878[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877))));
    }
    dist281867 = (((r81878[0] * r81878[0]) + (r81878[1] * r81878[1])) + (r81878[2] * r81878[2]));
    for (int loop_var81891 = 0; loop_var81891 < 24; loop_var81891++) {
      int mapped_idx;
      mapped_idx = loop_var81891;
      int al_index_name_symbol;
      al_index_name_symbol = ddist33107dinitial_mpg81805[loop_var81891];
      ddist33107dinitial_der81868[mapped_idx] = ((dist281867 * (0.5 * (sqrtl((1.0 / dist281867)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866)))) + (sqrtl(dist281867) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866)));
    }
    dist381869 = (dist281867 * sqrtl(dist281867));
    for (int slice_idx = 0; slice_idx < (((i81875 * 3) + 3) - (i81875 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81875 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81849[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81853[((24 * ((i81875 * 3) + slice_idx)) + mapped_idx)] = ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877) * -2.8247609439046209905e-07) * dist381869) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81806, ddist33107dinitial_der81868) * (r81878[(0 + slice_idx)] * -2.8247609439046209905e-07))) / (dist381869 * dist381869));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81875 * 3) + 3) - (i81875 * 3)); slice_idx++) {
      sat_acc81854[((i81875 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r81878[(0 + slice_idx)]) / dist381869);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81817[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81855[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81818, dcentral_acc2846dinitial_der81855) + ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81826, dr3376dinitial_der81877) * sat_gms81793[i81875]) * dist381869) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81806, ddist33107dinitial_der81868) * (r81878[(0 + slice_idx)] * sat_gms81793[i81875]))) / (dist381869 * dist381869)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81856[(0 + slice_idx)] = (central_acc81856[(0 + slice_idx)] + ((sat_gms81793[i81875] * r81878[(0 + slice_idx)]) / dist381869));
    }
  }
  for (int i81901 = 0; i81901 < 4; i81901++) {
    long double r81904[3] = { 0.0 };
    long double dr3708dinitial_der81903[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3708dinitial_mpg81827[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3708dinitial_der81903[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81904[(0 + slice_idx)] = (perturb_pos81792[((i81901 * 3) + slice_idx)] - central_pos81790[(0 + slice_idx)]);
    }
    for (int loop_var81911 = 0; loop_var81911 < 24; loop_var81911++) {
      int mapped_idx;
      mapped_idx = loop_var81911;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81813[loop_var81911];
      ddist23094dinitial_der81866[mapped_idx] = ((((r81904[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903))) + ((r81904[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)))) + ((r81904[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903))));
    }
    dist281867 = (((r81904[0] * r81904[0]) + (r81904[1] * r81904[1])) + (r81904[2] * r81904[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81817[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81855[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81818, dcentral_acc2846dinitial_der81855) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903) * perturb_gms81791[i81901]) * dist281867) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866) * (r81904[(0 + slice_idx)] * perturb_gms81791[i81901]))) / (dist281867 * dist281867)) * sqrtl(dist281867)) - ((0.5 * (sqrtl((1.0 / dist281867)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866))) * ((r81904[(0 + slice_idx)] * perturb_gms81791[i81901]) / dist281867))) / (sqrtl(dist281867) * sqrtl(dist281867))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81856[(0 + slice_idx)] = (central_acc81856[(0 + slice_idx)] + (((perturb_gms81791[i81901] * r81904[(0 + slice_idx)]) / dist281867) / sqrtl(dist281867)));
    }
    for (int j81917 = 0; j81917 < 4; j81917++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3708dinitial_mpg81827[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr3708dinitial_der81903[((24 * (0 + slice_idx)) + mapped_idx)] = -get_dfdx_cell(((j81917 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81904[(0 + slice_idx)] = (perturb_pos81792[((i81901 * 3) + slice_idx)] - (state81859[((j81917 * 6) + slice_idx)] + central_pos81790[(0 + slice_idx)]));
      }
      for (int loop_var81926 = 0; loop_var81926 < 24; loop_var81926++) {
        int mapped_idx;
        mapped_idx = loop_var81926;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81813[loop_var81926];
        ddist23094dinitial_der81866[mapped_idx] = ((((r81904[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903))) + ((r81904[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)))) + ((r81904[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903)) + (r81904[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903))));
      }
      dist281867 = (((r81904[0] * r81904[0]) + (r81904[1] * r81904[1])) + (r81904[2] * r81904[2]));
      for (int slice_idx = 0; slice_idx < (((j81917 * 3) + 3) - (j81917 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81917 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81849[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81853[((24 * ((j81917 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81917 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81850, dsat_acc2828dinitial_der81853) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81828, dr3708dinitial_der81903) * perturb_gms81791[i81901]) * dist281867) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866) * (r81904[(0 + slice_idx)] * perturb_gms81791[i81901]))) / (dist281867 * dist281867)) * sqrtl(dist281867)) - ((0.5 * (sqrtl((1.0 / dist281867)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866))) * ((r81904[(0 + slice_idx)] * perturb_gms81791[i81901]) / dist281867))) / (sqrtl(dist281867) * sqrtl(dist281867))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81917 * 3) + 3) - (j81917 * 3)); slice_idx++) {
        sat_acc81854[((j81917 * 3) + slice_idx)] = (sat_acc81854[((j81917 * 3) + slice_idx)] + (((perturb_gms81791[i81901] * r81904[(0 + slice_idx)]) / dist281867) / sqrtl(dist281867)));
      }
    }
  }
  for (int i81932 = 1; i81932 < 4; i81932++) {
    long double r81935[3] = { 0.0 };
    long double dr4219dinitial_der81934[72] = { 0.0 };
    for (int j81936 = 0; j81936 < i81932; j81936++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr4219dinitial_mpg81833[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr4219dinitial_der81934[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81936 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858) - get_dfdx_cell(((i81932 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81935[(0 + slice_idx)] = (state81859[((j81936 * 6) + slice_idx)] - state81859[((i81932 * 6) + slice_idx)]);
      }
      for (int loop_var81945 = 0; loop_var81945 < 24; loop_var81945++) {
        int mapped_idx;
        mapped_idx = loop_var81945;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81813[loop_var81945];
        ddist23094dinitial_der81866[mapped_idx] = ((((r81935[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934)) + (r81935[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934))) + ((r81935[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934)) + (r81935[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934)))) + ((r81935[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934)) + (r81935[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934))));
      }
      dist281867 = (((r81935[0] * r81935[0]) + (r81935[1] * r81935[1])) + (r81935[2] * r81935[2]));
      for (int slice_idx = 0; slice_idx < (((i81932 * 3) + 3) - (i81932 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((i81932 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81849[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81853[((24 * ((i81932 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81932 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81850, dsat_acc2828dinitial_der81853) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934) * sat_gms81793[j81936]) * dist281867) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866) * (r81935[(0 + slice_idx)] * sat_gms81793[j81936]))) / (dist281867 * dist281867)) * sqrtl(dist281867)) - ((0.5 * (sqrtl((1.0 / dist281867)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866))) * ((r81935[(0 + slice_idx)] * sat_gms81793[j81936]) / dist281867))) / (sqrtl(dist281867) * sqrtl(dist281867))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i81932 * 3) + 3) - (i81932 * 3)); slice_idx++) {
        sat_acc81854[((i81932 * 3) + slice_idx)] = (sat_acc81854[((i81932 * 3) + slice_idx)] + (((sat_gms81793[j81936] * r81935[(0 + slice_idx)]) / dist281867) / sqrtl(dist281867)));
      }
      for (int slice_idx = 0; slice_idx < (((j81936 * 3) + 3) - (j81936 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81936 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81849[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81853[((24 * ((j81936 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81936 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81850, dsat_acc2828dinitial_der81853) - (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81834, dr4219dinitial_der81934) * sat_gms81793[i81932]) * dist281867) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866) * (r81935[(0 + slice_idx)] * sat_gms81793[i81932]))) / (dist281867 * dist281867)) * sqrtl(dist281867)) - ((0.5 * (sqrtl((1.0 / dist281867)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81814, ddist23094dinitial_der81866))) * ((r81935[(0 + slice_idx)] * sat_gms81793[i81932]) / dist281867))) / (sqrtl(dist281867) * sqrtl(dist281867))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81936 * 3) + 3) - (j81936 * 3)); slice_idx++) {
        sat_acc81854[((j81936 * 3) + slice_idx)] = (sat_acc81854[((j81936 * 3) + slice_idx)] - (((sat_gms81793[i81932] * r81935[(0 + slice_idx)]) / dist281867) / sqrtl(dist281867)));
      }
    }
  }
  long double grav_acc_local81956[3] = { 0.0 };
  long double dgrav_acc_local4544dinitial_der81955[72] = { 0.0 };
  long double grav_acc81958[3] = { 0.0 };
  long double dgrav_acc4569dinitial_der81957[72] = { 0.0 };
  long double acc81960[3] = { 0.0 };
  long double dacc4588dinitial_der81959[72] = { 0.0 };
  long double rot81961[9] = { 0.0 };
  long double v_081963[7] = { 0.0 };
  long double dv_04620dinitial_der81962[168] = { 0.0 };
  long double dv_0_dx81965[7] = { 0.0 };
  long double ddv_0_dx4638dinitial_der81964[168] = { 0.0 };
  long double dv_0_dy81967[7] = { 0.0 };
  long double ddv_0_dy4660dinitial_der81966[168] = { 0.0 };
  long double dv_0_dz81969[7] = { 0.0 };
  long double ddv_0_dz4682dinitial_der81968[168] = { 0.0 };
  {
    long double jupiter_rotation_matrix81978[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix81978, t81788);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot81961[(0 + slice_idx)] = jupiter_rotation_matrix81978[slice_idx];
    }
  }
  for (int i81979 = 0; i81979 < 4; i81979++) {
    long double x81982;
    long double dx4794dinitial_der81981[24] = { 0.0 };
    for (int loop_var81986 = 0; loop_var81986 < 24; loop_var81986++) {
      int mapped_idx;
      mapped_idx = loop_var81986;
      int al_index_name_symbol;
      al_index_name_symbol = dx4794dinitial_mpg81801[loop_var81986];
      dx4794dinitial_der81981[mapped_idx] = (get_dfdx_cell(((i81979 * 6) + 0), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858) / 0.0004778945025452157572e0);
    }
    x81982 = (state81859[((i81979 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y81989;
    long double dy4843dinitial_der81988[24] = { 0.0 };
    for (int loop_var81993 = 0; loop_var81993 < 24; loop_var81993++) {
      int mapped_idx;
      mapped_idx = loop_var81993;
      int al_index_name_symbol;
      al_index_name_symbol = dy4843dinitial_mpg81829[loop_var81993];
      dy4843dinitial_der81988[mapped_idx] = (get_dfdx_cell(((i81979 * 6) + 1), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858) / 0.0004778945025452157572e0);
    }
    y81989 = (state81859[((i81979 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z81996;
    long double dz4892dinitial_der81995[24] = { 0.0 };
    for (int loop_var82000 = 0; loop_var82000 < 24; loop_var82000++) {
      int mapped_idx;
      mapped_idx = loop_var82000;
      int al_index_name_symbol;
      al_index_name_symbol = dz4892dinitial_mpg81823[loop_var82000];
      dz4892dinitial_der81995[mapped_idx] = (get_dfdx_cell(((i81979 * 6) + 2), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81800, dstate2924dinitial_der81858) / 0.0004778945025452157572e0);
    }
    z81996 = (state81859[((i81979 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82003;
    long double d_x4942dinitial_der82002[24] = { 0.0 };
    for (int loop_var82007 = 0; loop_var82007 < 24; loop_var82007++) {
      int mapped_idx;
      mapped_idx = loop_var82007;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4942dinitial_mpg81835[loop_var82007];
      d_x4942dinitial_der82002[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81802, dx4794dinitial_der81981) * rot81961[0]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81830, dy4843dinitial_der81988) * rot81961[3])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81824, dz4892dinitial_der81995) * rot81961[6]));
    }
    _x82003 = (((rot81961[0] * x81982) + (rot81961[3] * y81989)) + (rot81961[6] * z81996));
    long double _y82010;
    long double d_y4993dinitial_der82009[24] = { 0.0 };
    for (int loop_var82014 = 0; loop_var82014 < 24; loop_var82014++) {
      int mapped_idx;
      mapped_idx = loop_var82014;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4993dinitial_mpg81803[loop_var82014];
      d_y4993dinitial_der82009[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81802, dx4794dinitial_der81981) * rot81961[1]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81830, dy4843dinitial_der81988) * rot81961[4])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81824, dz4892dinitial_der81995) * rot81961[7]));
    }
    _y82010 = (((rot81961[1] * x81982) + (rot81961[4] * y81989)) + (rot81961[7] * z81996));
    long double _z82017;
    long double d_z5044dinitial_der82016[24] = { 0.0 };
    for (int loop_var82021 = 0; loop_var82021 < 24; loop_var82021++) {
      int mapped_idx;
      mapped_idx = loop_var82021;
      int al_index_name_symbol;
      al_index_name_symbol = d_z5044dinitial_mpg81807[loop_var82021];
      d_z5044dinitial_der82016[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81802, dx4794dinitial_der81981) * rot81961[2]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81830, dy4843dinitial_der81988) * rot81961[5])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81824, dz4892dinitial_der81995) * rot81961[8]));
    }
    _z82017 = (((rot81961[2] * x81982) + (rot81961[5] * y81989)) + (rot81961[8] * z81996));
    long double _r282024;
    long double d_r25100dinitial_der82023[24] = { 0.0 };
    for (int loop_var82028 = 0; loop_var82028 < 24; loop_var82028++) {
      int mapped_idx;
      mapped_idx = loop_var82028;
      int al_index_name_symbol;
      al_index_name_symbol = d_r25100dinitial_mpg81819[loop_var82028];
      d_r25100dinitial_der82023[mapped_idx] = ((((_x82003 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002)) + (_x82003 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002))) + ((_y82010 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009)) + (_y82010 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009)))) + ((_z82017 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016))));
    }
    _r282024 = (((_x82003 * _x82003) + (_y82010 * _y82010)) + (_z82017 * _z82017));
    long double _r82031;
    long double d_r5143dinitial_der82030[24] = { 0.0 };
    for (int loop_var82035 = 0; loop_var82035 < 24; loop_var82035++) {
      int mapped_idx;
      mapped_idx = loop_var82035;
      int al_index_name_symbol;
      al_index_name_symbol = d_r5143dinitial_mpg81839[loop_var82035];
      d_r5143dinitial_der82030[mapped_idx] = (0.5 * (sqrtl((1.0 / _r282024)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023)));
    }
    _r82031 = sqrtl(_r282024);
    long double _r382038;
    long double d_r35167dinitial_der82037[24] = { 0.0 };
    for (int loop_var82042 = 0; loop_var82042 < 24; loop_var82042++) {
      int mapped_idx;
      mapped_idx = loop_var82042;
      int al_index_name_symbol;
      al_index_name_symbol = d_r35167dinitial_mpg81821[loop_var82042];
      d_r35167dinitial_der82037[mapped_idx] = ((_r82031 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023)) + (_r282024 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81840, d_r5143dinitial_der82030)));
    }
    _r382038 = (_r82031 * _r282024);
    long double _r482045;
    long double d_r45191dinitial_der82044[24] = { 0.0 };
    for (int loop_var82049 = 0; loop_var82049 < 24; loop_var82049++) {
      int mapped_idx;
      mapped_idx = loop_var82049;
      int al_index_name_symbol;
      al_index_name_symbol = d_r45191dinitial_mpg81847[loop_var82049];
      d_r45191dinitial_der82044[mapped_idx] = ((_r282024 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023)) + (_r282024 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023)));
    }
    _r482045 = (_r282024 * _r282024);
    long double _r582052;
    long double d_r55216dinitial_der82051[24] = { 0.0 };
    for (int loop_var82056 = 0; loop_var82056 < 24; loop_var82056++) {
      int mapped_idx;
      mapped_idx = loop_var82056;
      int al_index_name_symbol;
      al_index_name_symbol = d_r55216dinitial_mpg81795[loop_var82056];
      d_r55216dinitial_der82051[mapped_idx] = ((_r482045 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81840, d_r5143dinitial_der82030)) + (_r82031 * get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044)));
    }
    _r582052 = (_r482045 * _r82031);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81845[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81962[((24 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r82031 * _r82031)) * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81840, d_r5143dinitial_der82030))));
        }
      }
    }
    v_081963[0] = (1.0e0 / _r82031);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81845[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81962[((24 * 1) + mapped_idx)] = ((((get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * 1.7320508075688772936e0) * _r382038) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81822, d_r35167dinitial_der82037) * (_z82017 * 1.7320508075688772936e0))) / (_r382038 * _r382038));
        }
      }
    }
    v_081963[1] = ((1.7320508075688772936e0 * _z82017) / _r382038);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81809[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81964[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002) * _r382038) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81822, d_r35167dinitial_der82037) * -_x82003)) / (_r382038 * _r382038));
        }
      }
    }
    dv_0_dx81965[0] = (-_x82003 / _r382038);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81841[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81966[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009) * _r382038) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81822, d_r35167dinitial_der82037) * -_y82010)) / (_r382038 * _r382038));
        }
      }
    }
    dv_0_dy81967[0] = (-_y82010 / _r382038);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81815[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81968[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * _r382038) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81822, d_r35167dinitial_der82037) * -_z82017)) / (_r382038 * _r382038));
        }
      }
    }
    dv_0_dz81969[0] = (-_z82017 / _r382038);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81809[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81964[((24 * 1) + mapped_idx)] = ((((((_z82017 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002)) + (_x82003 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * -5.1961524227066318805e0))) * _r582052) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81796, d_r55216dinitial_der82051) * ((_z82017 * -5.1961524227066318805e0) * _x82003))) / (_r582052 * _r582052));
        }
      }
    }
    dv_0_dx81965[1] = (((-5.1961524227066318805e0 * _z82017) * _x82003) / _r582052);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81841[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81966[((24 * 1) + mapped_idx)] = ((((((_z82017 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009)) + (_y82010 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * -5.1961524227066318805e0))) * _r582052) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81796, d_r55216dinitial_der82051) * ((_z82017 * -5.1961524227066318805e0) * _y82010))) / (_r582052 * _r582052));
        }
      }
    }
    dv_0_dy81967[1] = (((-5.1961524227066318805e0 * _z82017) * _y82010) / _r582052);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81815[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81968[((24 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r382038 * _r382038)) * get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81822, d_r35167dinitial_der82037)))) - ((((((_z82017 * 3.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * 3.0e0))) * _r582052) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81796, d_r55216dinitial_der82051) * ((_z82017 * 3.0e0) * _z82017))) / (_r582052 * _r582052))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz81969[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382038) - (((3.0e0 * _z82017) * _z82017) / _r582052)));
    for (int n82090 = 2; n82090 < 7; n82090++) {
      long double coef182092;
      coef182092 = sqrtl(((((2.0e0 * ((long double) n82090)) - 1.0e0) * ((long double) ((2 * n82090) + 1))) / ((long double) (n82090 * n82090))));
      long double coef282096;
      coef282096 = sqrtl(((((((long double) n82090) - 1.0e0) * ((long double) (n82090 - 1))) * ((long double) ((2 * n82090) + 1))) / ((long double) ((n82090 * n82090) * ((2 * n82090) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82090));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_04620dinitial_mpg81845[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dv_04620dinitial_der81962[((24 * n82090) + mapped_idx)] = (((((((v_081963[(n82090 - 1)] * coef182092) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * (get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * coef182092))) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * ((v_081963[(n82090 - 1)] * coef182092) * _z82017))) / (_r282024 * _r282024)) - ((((get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * coef282096) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * (v_081963[(n82090 - 2)] * coef282096))) / (_r282024 * _r282024)));
          }
        }
      }
      v_081963[n82090] = ((((coef182092 * v_081963[(n82090 - 1)]) * _z82017) / _r282024) - ((coef282096 * v_081963[(n82090 - 2)]) / _r282024));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82090));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx4638dinitial_mpg81809[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dx4638dinitial_der81964[((24 * n82090) + mapped_idx)] = ((((_z82017 * coef182092) * ((((get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81810, ddv_0_dx4638dinitial_der81964) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * dv_0_dx81965[(n82090 - 1)])) / (_r282024 * _r282024)) - ((((((v_081963[(n82090 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002)) + (_x82003 * (get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((v_081963[(n82090 - 1)] * 2.0e0) * _x82003))) / (_r482045 * _r482045)))) + (((dv_0_dx81965[(n82090 - 1)] / _r282024) - (((v_081963[(n82090 - 1)] * 2.0e0) * _x82003) / _r482045)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * coef182092))) - (((((get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81810, ddv_0_dx4638dinitial_der81964) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * dv_0_dx81965[(n82090 - 2)])) / (_r282024 * _r282024)) - ((((((v_081963[(n82090 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81836, d_x4942dinitial_der82002)) + (_x82003 * (get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((v_081963[(n82090 - 2)] * 2.0e0) * _x82003))) / (_r482045 * _r482045))) * coef282096));
          }
        }
      }
      dv_0_dx81965[n82090] = (((coef182092 * _z82017) * ((dv_0_dx81965[(n82090 - 1)] / _r282024) - (((v_081963[(n82090 - 1)] * 2.0e0) * _x82003) / _r482045))) - (coef282096 * ((dv_0_dx81965[(n82090 - 2)] / _r282024) - (((v_081963[(n82090 - 2)] * 2.0e0) * _x82003) / _r482045))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82090));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy4660dinitial_mpg81841[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dy4660dinitial_der81966[((24 * n82090) + mapped_idx)] = ((((_z82017 * coef182092) * ((((get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81842, ddv_0_dy4660dinitial_der81966) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * dv_0_dy81967[(n82090 - 1)])) / (_r282024 * _r282024)) - ((((((v_081963[(n82090 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009)) + (_y82010 * (get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((v_081963[(n82090 - 1)] * 2.0e0) * _y82010))) / (_r482045 * _r482045)))) + (((dv_0_dy81967[(n82090 - 1)] / _r282024) - (((v_081963[(n82090 - 1)] * 2.0e0) * _y82010) / _r482045)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * coef182092))) - (((((get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81842, ddv_0_dy4660dinitial_der81966) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * dv_0_dy81967[(n82090 - 2)])) / (_r282024 * _r282024)) - ((((((v_081963[(n82090 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81804, d_y4993dinitial_der82009)) + (_y82010 * (get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((v_081963[(n82090 - 2)] * 2.0e0) * _y82010))) / (_r482045 * _r482045))) * coef282096));
          }
        }
      }
      dv_0_dy81967[n82090] = (((coef182092 * _z82017) * ((dv_0_dy81967[(n82090 - 1)] / _r282024) - (((v_081963[(n82090 - 1)] * 2.0e0) * _y82010) / _r482045))) - (coef282096 * ((dv_0_dy81967[(n82090 - 2)] / _r282024) - (((v_081963[(n82090 - 2)] * 2.0e0) * _y82010) / _r482045))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82090));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz4682dinitial_mpg81815[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dz4682dinitial_der81968[((24 * n82090) + mapped_idx)] = ((((((((dv_0_dz81969[(n82090 - 1)] * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81816, ddv_0_dz4682dinitial_der81968))) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * (dv_0_dz81969[(n82090 - 1)] * _z82017))) / (_r282024 * _r282024)) + ((v_081963[(n82090 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r282024 * _r282024)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023)))) - ((((((_z82017 * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((_z82017 * 2.0e0) * _z82017))) / (_r482045 * _r482045)))) + (((1.0e0 / _r282024) - (((_z82017 * 2.0e0) * _z82017) / _r482045)) * get_dfdx_cell((n82090 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962)))) * coef182092) - (((((get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81816, ddv_0_dz4682dinitial_der81968) * _r282024) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81820, d_r25100dinitial_der82023) * dv_0_dz81969[(n82090 - 2)])) / (_r282024 * _r282024)) - ((((((v_081963[(n82090 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81808, d_z5044dinitial_der82016)) + (_z82017 * (get_dfdx_cell((n82090 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81846, dv_04620dinitial_der81962) * 2.0e0))) * _r482045) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81848, d_r45191dinitial_der82044) * ((v_081963[(n82090 - 2)] * 2.0e0) * _z82017))) / (_r482045 * _r482045))) * coef282096));
          }
        }
      }
      dv_0_dz81969[n82090] = ((coef182092 * (((dv_0_dz81969[(n82090 - 1)] * _z82017) / _r282024) + (v_081963[(n82090 - 1)] * ((1.0e0 / _r282024) - (((2.0e0 * _z82017) * _z82017) / _r482045))))) - (coef282096 * ((dv_0_dz81969[(n82090 - 2)] / _r282024) - (((v_081963[(n82090 - 2)] * 2.0e0) * _z82017) / _r482045))));
    }
    long double dpotential_dx82117;
    long double ddpotential_dx6387dinitial_der82116[24] = { 0.0 };
    for (int loop_var82121 = 0; loop_var82121 < 24; loop_var82121++) {
      int mapped_idx;
      mapped_idx = loop_var82121;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx6387dinitial_mpg81811[loop_var82121];
      ddpotential_dx6387dinitial_der82116[mapped_idx] = 0.0;
    }
    dpotential_dx82117 = ((long double) 0);
    long double dpotential_dy82123;
    long double ddpotential_dy6414dinitial_der82122[24] = { 0.0 };
    for (int loop_var82127 = 0; loop_var82127 < 24; loop_var82127++) {
      int mapped_idx;
      mapped_idx = loop_var82127;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy6414dinitial_mpg81851[loop_var82127];
      ddpotential_dy6414dinitial_der82122[mapped_idx] = 0.0;
    }
    dpotential_dy82123 = ((long double) 0);
    long double dpotential_dz82129;
    long double ddpotential_dz6441dinitial_der82128[24] = { 0.0 };
    for (int loop_var82133 = 0; loop_var82133 < 24; loop_var82133++) {
      int mapped_idx;
      mapped_idx = loop_var82133;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz6441dinitial_mpg81843[loop_var82133];
      ddpotential_dz6441dinitial_der82128[mapped_idx] = 0.0;
    }
    dpotential_dz82129 = ((long double) 0);
    for (int n82134 = 2; n82134 < 7; n82134++) {
      for (int loop_var82139 = 0; loop_var82139 < 24; loop_var82139++) {
        int mapped_idx;
        mapped_idx = loop_var82139;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx6387dinitial_mpg81811[loop_var82139];
        ddpotential_dx6387dinitial_der82116[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81812, ddpotential_dx6387dinitial_der82116) + (get_dfdx_cell(n82134, 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81810, ddv_0_dx4638dinitial_der81964) * central_grav81722[n82134]));
      }
      dpotential_dx82117 = (dpotential_dx82117 + (dv_0_dx81965[n82134] * central_grav81722[n82134]));
      for (int loop_var82144 = 0; loop_var82144 < 24; loop_var82144++) {
        int mapped_idx;
        mapped_idx = loop_var82144;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy6414dinitial_mpg81851[loop_var82144];
        ddpotential_dy6414dinitial_der82122[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81852, ddpotential_dy6414dinitial_der82122) + (get_dfdx_cell(n82134, 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81842, ddv_0_dy4660dinitial_der81966) * central_grav81722[n82134]));
      }
      dpotential_dy82123 = (dpotential_dy82123 + (dv_0_dy81967[n82134] * central_grav81722[n82134]));
      for (int loop_var82149 = 0; loop_var82149 < 24; loop_var82149++) {
        int mapped_idx;
        mapped_idx = loop_var82149;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz6441dinitial_mpg81843[loop_var82149];
        ddpotential_dz6441dinitial_der82128[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81844, ddpotential_dz6441dinitial_der82128) + (get_dfdx_cell(n82134, 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81816, ddv_0_dz4682dinitial_der81968) * central_grav81722[n82134]));
      }
      dpotential_dz82129 = (dpotential_dz82129 + (dv_0_dz81969[n82134] * central_grav81722[n82134]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81831[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81955[((24 * 0) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81812, ddpotential_dx6387dinitial_der82116);
        }
      }
    }
    grav_acc_local81956[0] = dpotential_dx82117;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81831[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81955[((24 * 1) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81852, ddpotential_dy6414dinitial_der82122);
        }
      }
    }
    grav_acc_local81956[1] = dpotential_dy82123;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 2));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81831[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81955[((24 * 2) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81844, ddpotential_dz6441dinitial_der82128);
        }
      }
    }
    grav_acc_local81956[2] = dpotential_dz82129;
    for (int k82163 = 0; k82163 < 3; k82163++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * k82163));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc4569dinitial_mpg81837[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dgrav_acc4569dinitial_der81957[((24 * k82163) + mapped_idx)] = ((((get_dfdx_cell(0, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81832, dgrav_acc_local4544dinitial_der81955) * rot81961[((3 * k82163) + 0)]) + (get_dfdx_cell(1, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81832, dgrav_acc_local4544dinitial_der81955) * rot81961[((3 * k82163) + 1)])) + (get_dfdx_cell(2, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81832, dgrav_acc_local4544dinitial_der81955) * rot81961[((3 * k82163) + 2)])) / 2.2838315556293922983e-07);
          }
        }
      }
      grav_acc81958[k82163] = ((((rot81961[((3 * k82163) + 0)] * grav_acc_local81956[0]) + (rot81961[((3 * k82163) + 1)] * grav_acc_local81956[1])) + (rot81961[((3 * k82163) + 2)] * grav_acc_local81956[2])) / 2.2838315556293922983e-07);
    }
    for (int slice_idx = 0; slice_idx < (((i81979 * 3) + 3) - (i81979 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81979 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81849[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81853[((24 * ((i81979 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81979 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81850, dsat_acc2828dinitial_der81853) + (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81838, dgrav_acc4569dinitial_der81957) * 2.8247609439046209905e-07));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81979 * 3) + 3) - (i81979 * 3)); slice_idx++) {
      sat_acc81854[((i81979 * 3) + slice_idx)] = (sat_acc81854[((i81979 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc81958[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81817[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81855[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81818, dcentral_acc2846dinitial_der81855) - (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81838, dgrav_acc4569dinitial_der81957) * sat_gms81793[i81979]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81856[(0 + slice_idx)] = (central_acc81856[(0 + slice_idx)] - (sat_gms81793[i81979] * grav_acc81958[(0 + slice_idx)]));
    }
  }
  for (int i82177 = 0; i82177 < 4; i82177++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc4588dinitial_mpg81797[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dacc4588dinitial_der81959[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82177 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81850, dsat_acc2828dinitial_der81853) - get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81818, dcentral_acc2846dinitial_der81855));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc81960[(0 + slice_idx)] = (sat_acc81854[((i82177 * 3) + slice_idx)] - central_acc81856[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82177 * 6) + 3) - (i82177 * 6)); slice_idx++) {
        jupsatsystem81787[((i82177 * 6) + slice_idx)] = state81859[(((i82177 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82177 * 6) + 6) - ((i82177 * 6) + 3)); slice_idx++) {
        jupsatsystem81787[(((i82177 * 6) + 3) + slice_idx)] = acc81960[(0 + slice_idx)];
      }
    }
    for (int j82189 = 0; j82189 < 3; j82189++) {
      {
        for (int slice_idx = 0; slice_idx < ((24 + (((((i82177 * 6) + j82189) + 1) * 4) * 6)) - (24 + ((((i82177 * 6) + j82189) * 4) * 6))); slice_idx++) {
          jupsatsystem81787[((24 + ((((i82177 * 6) + j82189) * 4) * 6)) + slice_idx)] = state_and_derivatives81789[((24 + (((((i82177 * 6) + j82189) + 3) * 4) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((24 + (((((i82177 * 6) + j82189) + 4) * 4) * 6)) - (24 + (((((i82177 * 6) + j82189) + 3) * 4) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (24 + (((((i82177 * 6) + j82189) + 3) * 4) * 6)));
        jupsatsystem81787[func_slice_idx] = get_dfdx_cell(j82189, 24, al_index_name_symbol, 24, dacc4588dinitial_inv_mpg81798, dacc4588dinitial_der81959);
      }
    }
  }
  return 0;
}

int jupsatsystem_noderiv(long double *restrict jupsatsystem_noderiv82196, long double t82197, long double *restrict state82198, long double *restrict central_pos82199, long double *restrict perturb_gms82200, long double *restrict perturb_pos82201, long double *restrict sat_gms82202) {
  long double sat_acc82204[12] = { 0.0 };
  long double central_acc82205[3] = { 0.0 };
  long double dist282206;
  long double dist382207;
  for (int i82208 = 0; i82208 < 4; i82208++) {
    long double r82210[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82210[(0 + slice_idx)] = state82198[((i82208 * 6) + slice_idx)];
      }
    }
    dist282206 = (((r82210[0] * r82210[0]) + (r82210[1] * r82210[1])) + (r82210[2] * r82210[2]));
    dist382207 = (dist282206 * sqrtl(dist282206));
    {
      for (int slice_idx = 0; slice_idx < (((i82208 * 3) + 3) - (i82208 * 3)); slice_idx++) {
        sat_acc82204[((i82208 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r82210[(0 + slice_idx)]) / dist382207);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82205[(0 + slice_idx)] = (central_acc82205[(0 + slice_idx)] + ((sat_gms82202[i82208] * r82210[(0 + slice_idx)]) / dist382207));
      }
    }
  }
  for (int i82226 = 0; i82226 < 4; i82226++) {
    long double r82228[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82228[(0 + slice_idx)] = (perturb_pos82201[((i82226 * 3) + slice_idx)] - central_pos82199[(0 + slice_idx)]);
      }
    }
    dist282206 = (((r82228[0] * r82228[0]) + (r82228[1] * r82228[1])) + (r82228[2] * r82228[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82205[(0 + slice_idx)] = (central_acc82205[(0 + slice_idx)] + (((perturb_gms82200[i82226] * r82228[(0 + slice_idx)]) / dist282206) / sqrtl(dist282206)));
      }
    }
    for (int j82238 = 0; j82238 < 4; j82238++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82228[(0 + slice_idx)] = (perturb_pos82201[((i82226 * 3) + slice_idx)] - (state82198[((j82238 * 6) + slice_idx)] + central_pos82199[(0 + slice_idx)]));
        }
      }
      dist282206 = (((r82228[0] * r82228[0]) + (r82228[1] * r82228[1])) + (r82228[2] * r82228[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j82238 * 3) + 3) - (j82238 * 3)); slice_idx++) {
          sat_acc82204[((j82238 * 3) + slice_idx)] = (sat_acc82204[((j82238 * 3) + slice_idx)] + (((perturb_gms82200[i82226] * r82228[(0 + slice_idx)]) / dist282206) / sqrtl(dist282206)));
        }
      }
    }
  }
  for (int i82249 = 1; i82249 < 4; i82249++) {
    long double r82251[3] = { 0.0 };
    for (int j82252 = 0; j82252 < i82249; j82252++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82251[(0 + slice_idx)] = (state82198[((j82252 * 6) + slice_idx)] - state82198[((i82249 * 6) + slice_idx)]);
        }
      }
      dist282206 = (((r82251[0] * r82251[0]) + (r82251[1] * r82251[1])) + (r82251[2] * r82251[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i82249 * 3) + 3) - (i82249 * 3)); slice_idx++) {
          sat_acc82204[((i82249 * 3) + slice_idx)] = (sat_acc82204[((i82249 * 3) + slice_idx)] + (((sat_gms82202[j82252] * r82251[(0 + slice_idx)]) / dist282206) / sqrtl(dist282206)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j82252 * 3) + 3) - (j82252 * 3)); slice_idx++) {
          sat_acc82204[((j82252 * 3) + slice_idx)] = (sat_acc82204[((j82252 * 3) + slice_idx)] - (((sat_gms82202[i82249] * r82251[(0 + slice_idx)]) / dist282206) / sqrtl(dist282206)));
        }
      }
    }
  }
  long double grav_acc_local82266[3] = { 0.0 };
  long double grav_acc82267[3] = { 0.0 };
  long double acc82268[3] = { 0.0 };
  long double rot82269[9] = { 0.0 };
  long double v_082270[7] = { 0.0 };
  long double dv_0_dx82271[7] = { 0.0 };
  long double dv_0_dy82272[7] = { 0.0 };
  long double dv_0_dz82273[7] = { 0.0 };
  {
    long double jupiter_rotation_matrix82282[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix82282, t82197);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot82269[(0 + slice_idx)] = jupiter_rotation_matrix82282[slice_idx];
    }
  }
  for (int i82283 = 0; i82283 < 4; i82283++) {
    long double x82285;
    x82285 = (state82198[((i82283 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y82289;
    y82289 = (state82198[((i82283 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z82293;
    z82293 = (state82198[((i82283 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82297;
    _x82297 = (((rot82269[0] * x82285) + (rot82269[3] * y82289)) + (rot82269[6] * z82293));
    long double _y82301;
    _y82301 = (((rot82269[1] * x82285) + (rot82269[4] * y82289)) + (rot82269[7] * z82293));
    long double _z82305;
    _z82305 = (((rot82269[2] * x82285) + (rot82269[5] * y82289)) + (rot82269[8] * z82293));
    long double _r282309;
    _r282309 = (((_x82297 * _x82297) + (_y82301 * _y82301)) + (_z82305 * _z82305));
    long double r82313;
    r82313 = sqrtl(_r282309);
    long double _r382317;
    _r382317 = (r82313 * _r282309);
    long double _r482321;
    _r482321 = (_r282309 * _r282309);
    long double _r582325;
    _r582325 = (_r482321 * r82313);
    v_082270[0] = (1.0e0 / r82313);
    v_082270[1] = ((1.7320508075688772936e0 * _z82305) / _r382317);
    dv_0_dx82271[0] = (-_x82297 / _r382317);
    dv_0_dy82272[0] = (-_y82301 / _r382317);
    dv_0_dz82273[0] = (-_z82305 / _r382317);
    dv_0_dx82271[1] = (((-5.1961524227066318805e0 * _z82305) * _x82297) / _r582325);
    dv_0_dy82272[1] = (((-5.1961524227066318805e0 * _z82305) * _y82301) / _r582325);
    dv_0_dz82273[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382317) - (((3.0e0 * _z82305) * _z82305) / _r582325)));
    for (int n82353 = 2; n82353 < 7; n82353++) {
      long double coef182355;
      coef182355 = sqrtl(((((2.0e0 * ((long double) n82353)) - 1.0e0) * ((long double) ((2 * n82353) + 1))) / ((long double) (n82353 * n82353))));
      long double coef282359;
      coef282359 = sqrtl(((((((long double) n82353) - 1.0e0) * ((long double) (n82353 - 1))) * ((long double) ((2 * n82353) + 1))) / ((long double) ((n82353 * n82353) * ((2 * n82353) - 3)))));
      v_082270[n82353] = ((((coef182355 * v_082270[(n82353 - 1)]) * _z82305) / _r282309) - ((coef282359 * v_082270[(n82353 - 2)]) / _r282309));
      dv_0_dx82271[n82353] = (((coef182355 * _z82305) * ((dv_0_dx82271[(n82353 - 1)] / _r282309) - (((v_082270[(n82353 - 1)] * 2.0e0) * _x82297) / _r482321))) - (coef282359 * ((dv_0_dx82271[(n82353 - 2)] / _r282309) - (((v_082270[(n82353 - 2)] * 2.0e0) * _x82297) / _r482321))));
      dv_0_dy82272[n82353] = (((coef182355 * _z82305) * ((dv_0_dy82272[(n82353 - 1)] / _r282309) - (((v_082270[(n82353 - 1)] * 2.0e0) * _y82301) / _r482321))) - (coef282359 * ((dv_0_dy82272[(n82353 - 2)] / _r282309) - (((v_082270[(n82353 - 2)] * 2.0e0) * _y82301) / _r482321))));
      dv_0_dz82273[n82353] = ((coef182355 * (((dv_0_dz82273[(n82353 - 1)] * _z82305) / _r282309) + (v_082270[(n82353 - 1)] * ((1.0e0 / _r282309) - (((2.0e0 * _z82305) * _z82305) / _r482321))))) - (coef282359 * ((dv_0_dz82273[(n82353 - 2)] / _r282309) - (((v_082270[(n82353 - 2)] * 2.0e0) * _z82305) / _r482321))));
    }
    long double dpotential_dx82375;
    dpotential_dx82375 = ((long double) 0);
    long double dpotential_dy82379;
    dpotential_dy82379 = ((long double) 0);
    long double dpotential_dz82383;
    dpotential_dz82383 = ((long double) 0);
    for (int n82387 = 2; n82387 < 7; n82387++) {
      dpotential_dx82375 = (dpotential_dx82375 + (dv_0_dx82271[n82387] * central_grav81722[n82387]));
      dpotential_dy82379 = (dpotential_dy82379 + (dv_0_dy82272[n82387] * central_grav81722[n82387]));
      dpotential_dz82383 = (dpotential_dz82383 + (dv_0_dz82273[n82387] * central_grav81722[n82387]));
    }
    grav_acc_local82266[0] = dpotential_dx82375;
    grav_acc_local82266[1] = dpotential_dy82379;
    grav_acc_local82266[2] = dpotential_dz82383;
    for (int k82407 = 0; k82407 < 3; k82407++) {
      grav_acc82267[k82407] = ((((rot82269[((3 * k82407) + 0)] * grav_acc_local82266[0]) + (rot82269[((3 * k82407) + 1)] * grav_acc_local82266[1])) + (rot82269[((3 * k82407) + 2)] * grav_acc_local82266[2])) / 2.2838315556293922983e-07);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82283 * 3) + 3) - (i82283 * 3)); slice_idx++) {
        sat_acc82204[((i82283 * 3) + slice_idx)] = (sat_acc82204[((i82283 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc82267[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82205[(0 + slice_idx)] = (central_acc82205[(0 + slice_idx)] - (sat_gms82202[i82283] * grav_acc82267[(0 + slice_idx)]));
      }
    }
  }
  for (int i82418 = 0; i82418 < 4; i82418++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc82268[(0 + slice_idx)] = (sat_acc82204[((i82418 * 3) + slice_idx)] - central_acc82205[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82418 * 6) + 3) - (i82418 * 6)); slice_idx++) {
        jupsatsystem_noderiv82196[((i82418 * 6) + slice_idx)] = state82198[(((i82418 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82418 * 6) + 6) - ((i82418 * 6) + 3)); slice_idx++) {
        jupsatsystem_noderiv82196[(((i82418 * 6) + 3) + slice_idx)] = acc82268[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

