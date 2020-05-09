#include <math.h>

static inline long double get_dfdx_cell(int full_idx81823, int dx_mapped_size81824, int al_idx81825, int inv_mapping_period81826, int *restrict dx_idx_mappings81827, long double *restrict der_vec81828) {
  return ((al_idx81825 < inv_mapping_period81826) ? ((dx_idx_mappings81827[((inv_mapping_period81826 * full_idx81823) + al_idx81825)] >= 0) ? der_vec81828[((dx_mapped_size81824 * full_idx81823) + dx_idx_mappings81827[((inv_mapping_period81826 * full_idx81823) + al_idx81825)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx81823, int dx_mapped_size81824, int al_idx81825, int inv_mapping_period81826, int *restrict dx_idx_mappings81827, long double *restrict der_vec81828) {
  return ((al_idx81825 < inv_mapping_period81826) ? ((al_idx81825 == full_idx81823) ? 1.0 : ((dx_idx_mappings81827[((inv_mapping_period81826 * full_idx81823) + al_idx81825)] >= 0) ? der_vec81828[((dx_mapped_size81824 * full_idx81823) + dx_idx_mappings81827[((inv_mapping_period81826 * full_idx81823) + al_idx81825)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx81830, int inv_mapping_period81831, int *restrict dx_idx_mappings81832, long double *restrict der_vec81833) {
  return ((((al_idx81830 < inv_mapping_period81831) ? dx_idx_mappings81832[al_idx81830] : -1) >= 0) ? der_vec81833[((al_idx81830 < inv_mapping_period81831) ? dx_idx_mappings81832[al_idx81830] : -1)] : 0.0);
}

const long double pole_ra81834[2] = { 268.05659500000001572e0, -0.006498999999999999569e0 };
const long double pole_dec81835[2] = { 64.49530300000000693e0, 0.0024130000000000002142e0 };
const long double pm81836[2] = { 284.94999999999998863e0, 870.5359999999999445e0 };
const long double nut_prec_ra81837[5] = { 0.00011699999999999999788e0, 0.00093800000000000003167e0, 0.0014319999999999998962e0, 3.000000000000000076e-05, 0.0021500000000000000014e0 };
const long double nut_prec_dec81838[5] = { 5.0000000000000002396e-05, 0.00040400000000000000798e0, 0.0006170000000000000354e0, -1.29999999999999992e-05, 0.0009259999999999999568e0 };
const long double Jabcde_081839[5] = { 99.36071400000000153e0, 175.89536899999998809e0, 300.32316200000002482e0, 114.01230499999999779e0, 49.511251000000001454e0 };
const long double Jabcde_T81840[5] = { 4850.4045999999998457e0, 1191.9604999999999109e0, 262.54750000000001364e0, 6070.247599999999693e0, 64.29999999999999716e0 };
const long double central_grav81841[7] = { 0.0e0, 0.0e0, -0.0065724808672554692115e0, 0.0e0, 0.00019554099999999999462e0, 0.0e0, -9.4975767597683728686e-06 };

int jupiter_rotation_matrix(long double *restrict jupiter_rotation_matrix81842, long double t81843) {
  long double T81845;
  T81845 = (t81843 / 36525.0e0);
  long double alpha_081849;
  alpha_081849 = (268.05659500000001572e0 + (-0.006498999999999999569e0 * T81845));
  long double delta_081853;
  delta_081853 = (64.49530300000000693e0 + (0.0024130000000000002142e0 * T81845));
  long double W81857;
  W81857 = (((284.94999999999998863e0 + (870.5359999999999445e0 * t81843)) * 3.141592653589793116e0) / 180.0e0);
  for (int i81861 = 0; i81861 < 5; i81861++) {
    long double J81863;
    J81863 = (((Jabcde_081839[i81861] + (Jabcde_T81840[i81861] * T81845)) * 3.141592653589793116e0) / 180.0e0);
    alpha_081849 = (alpha_081849 + (nut_prec_ra81837[i81861] * sinl(J81863)));
    delta_081853 = (delta_081853 + (nut_prec_dec81838[i81861] * cosl(J81863)));
  }
  alpha_081849 = (alpha_081849 * 0.017453292519943295088e0);
  delta_081853 = (delta_081853 * 0.017453292519943295088e0);
  jupiter_rotation_matrix81842[0] = ((-sinl(alpha_081849) * cosl(W81857)) - ((cosl(alpha_081849) * sinl(delta_081853)) * sinl(W81857)));
  jupiter_rotation_matrix81842[1] = ((sinl(alpha_081849) * sinl(W81857)) - ((cosl(alpha_081849) * sinl(delta_081853)) * cosl(W81857)));
  jupiter_rotation_matrix81842[2] = (cosl(alpha_081849) * cosl(delta_081853));
  jupiter_rotation_matrix81842[3] = ((cosl(alpha_081849) * cosl(W81857)) - ((sinl(alpha_081849) * sinl(delta_081853)) * sinl(W81857)));
  jupiter_rotation_matrix81842[4] = ((-cosl(alpha_081849) * sinl(W81857)) - ((sinl(alpha_081849) * sinl(delta_081853)) * cosl(W81857)));
  jupiter_rotation_matrix81842[5] = (cosl(delta_081853) * sinl(alpha_081849));
  jupiter_rotation_matrix81842[6] = (cosl(delta_081853) * sinl(W81857));
  jupiter_rotation_matrix81842[7] = (cosl(delta_081853) * cosl(W81857));
  jupiter_rotation_matrix81842[8] = sinl(delta_081853);
  return 0;
}

int jupsatsystem(long double *restrict jupsatsystem81906, long double t81907, long double *restrict state_and_derivatives81908, long double *restrict central_pos81909, long double *restrict perturb_gms81910, long double *restrict perturb_pos81911, long double *restrict sat_gms81912) {
  static int d_r55216dinitial_mpg81914[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_inv_mpg81915[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_mpg81916[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_inv_mpg81917[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dstate2924dinitial_mpg81918[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2924dinitial_inv_mpg81919[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dx4794dinitial_mpg81920[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_inv_mpg81921[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_mpg81922[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_inv_mpg81923[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_mpg81924[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_inv_mpg81925[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_mpg81926[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_inv_mpg81927[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_mpg81928[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_inv_mpg81929[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_mpg81930[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_inv_mpg81931[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_mpg81932[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_inv_mpg81933[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_mpg81934[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_inv_mpg81935[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_mpg81936[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_inv_mpg81937[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_mpg81938[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_inv_mpg81939[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_mpg81940[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_inv_mpg81941[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_mpg81942[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_inv_mpg81943[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_mpg81944[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_inv_mpg81945[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_mpg81946[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_inv_mpg81947[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_mpg81948[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_inv_mpg81949[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_mpg81950[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_inv_mpg81951[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_mpg81952[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_inv_mpg81953[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_mpg81954[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_inv_mpg81955[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_mpg81956[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_inv_mpg81957[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_mpg81958[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_inv_mpg81959[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_mpg81960[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_inv_mpg81961[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_mpg81962[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_inv_mpg81963[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dv_04620dinitial_mpg81964[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_04620dinitial_inv_mpg81965[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int d_r45191dinitial_mpg81966[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_inv_mpg81967[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_mpg81968[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_inv_mpg81969[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_mpg81970[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_inv_mpg81971[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  long double sat_acc81973[12] = { 0.0 };
  long double dsat_acc2828dinitial_der81972[288] = { 0.0 };
  long double central_acc81975[3] = { 0.0 };
  long double dcentral_acc2846dinitial_der81974[72] = { 0.0 };
  long double state_derivatives_initial81976[576] = { 0.0 };
  long double state81978[24] = { 0.0 };
  long double dstate2924dinitial_der81977[576] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 576; slice_idx++) {
      state_derivatives_initial81976[(0 + slice_idx)] = state_and_derivatives81908[(24 + slice_idx)];
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
        al_index_name_symbol = dstate2924dinitial_mpg81918[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dstate2924dinitial_der81977[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    state81978[(0 + slice_idx)] = state_and_derivatives81908[(0 + slice_idx)];
  }
  long double dist281986;
  long double ddist23094dinitial_der81985[24] = { 0.0 };
  long double dist381988;
  long double ddist33107dinitial_der81987[24] = { 0.0 };
  for (int i81989 = 0; i81989 < 24; i81989++) {
    for (int j81991 = 0; j81991 < 24; j81991++) {
      {
        int mapped_idx;
        mapped_idx = ((j81991 < 24) ? dstate2924dinitial_inv_mpg81919[((i81989 * 24) + j81991)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2924dinitial_der81977[((i81989 * 24) + mapped_idx)] = state_derivatives_initial81976[(((i81989 * 4) * 6) + j81991)];
        }
      }
    }
  }
  for (int i81994 = 0; i81994 < 4; i81994++) {
    long double r81997[3] = { 0.0 };
    long double dr3376dinitial_der81996[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3376dinitial_mpg81944[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3376dinitial_der81996[((24 * (0 + slice_idx)) + mapped_idx)] = get_dfdx_cell(((i81994 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81997[(0 + slice_idx)] = state81978[((i81994 * 6) + slice_idx)];
    }
    for (int loop_var82005 = 0; loop_var82005 < 24; loop_var82005++) {
      int mapped_idx;
      mapped_idx = loop_var82005;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81932[loop_var82005];
      ddist23094dinitial_der81985[mapped_idx] = ((((r81997[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996)) + (r81997[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996))) + ((r81997[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996)) + (r81997[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996)))) + ((r81997[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996)) + (r81997[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996))));
    }
    dist281986 = (((r81997[0] * r81997[0]) + (r81997[1] * r81997[1])) + (r81997[2] * r81997[2]));
    for (int loop_var82010 = 0; loop_var82010 < 24; loop_var82010++) {
      int mapped_idx;
      mapped_idx = loop_var82010;
      int al_index_name_symbol;
      al_index_name_symbol = ddist33107dinitial_mpg81924[loop_var82010];
      ddist33107dinitial_der81987[mapped_idx] = ((dist281986 * (0.5 * (sqrtl((1.0 / dist281986)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985)))) + (sqrtl(dist281986) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985)));
    }
    dist381988 = (dist281986 * sqrtl(dist281986));
    for (int slice_idx = 0; slice_idx < (((i81994 * 3) + 3) - (i81994 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81994 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81968[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81972[((24 * ((i81994 * 3) + slice_idx)) + mapped_idx)] = ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996) * -2.8247609439046209905e-07) * dist381988) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81925, ddist33107dinitial_der81987) * (r81997[(0 + slice_idx)] * -2.8247609439046209905e-07))) / (dist381988 * dist381988));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81994 * 3) + 3) - (i81994 * 3)); slice_idx++) {
      sat_acc81973[((i81994 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r81997[(0 + slice_idx)]) / dist381988);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81936[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81974[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81937, dcentral_acc2846dinitial_der81974) + ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81945, dr3376dinitial_der81996) * sat_gms81912[i81994]) * dist381988) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81925, ddist33107dinitial_der81987) * (r81997[(0 + slice_idx)] * sat_gms81912[i81994]))) / (dist381988 * dist381988)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81975[(0 + slice_idx)] = (central_acc81975[(0 + slice_idx)] + ((sat_gms81912[i81994] * r81997[(0 + slice_idx)]) / dist381988));
    }
  }
  for (int i82020 = 0; i82020 < 4; i82020++) {
    long double r82023[3] = { 0.0 };
    long double dr3708dinitial_der82022[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3708dinitial_mpg81946[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3708dinitial_der82022[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r82023[(0 + slice_idx)] = (perturb_pos81911[((i82020 * 3) + slice_idx)] - central_pos81909[(0 + slice_idx)]);
    }
    for (int loop_var82030 = 0; loop_var82030 < 24; loop_var82030++) {
      int mapped_idx;
      mapped_idx = loop_var82030;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81932[loop_var82030];
      ddist23094dinitial_der81985[mapped_idx] = ((((r82023[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022))) + ((r82023[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)))) + ((r82023[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022))));
    }
    dist281986 = (((r82023[0] * r82023[0]) + (r82023[1] * r82023[1])) + (r82023[2] * r82023[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81936[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81974[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81937, dcentral_acc2846dinitial_der81974) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022) * perturb_gms81910[i82020]) * dist281986) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985) * (r82023[(0 + slice_idx)] * perturb_gms81910[i82020]))) / (dist281986 * dist281986)) * sqrtl(dist281986)) - ((0.5 * (sqrtl((1.0 / dist281986)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985))) * ((r82023[(0 + slice_idx)] * perturb_gms81910[i82020]) / dist281986))) / (sqrtl(dist281986) * sqrtl(dist281986))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81975[(0 + slice_idx)] = (central_acc81975[(0 + slice_idx)] + (((perturb_gms81910[i82020] * r82023[(0 + slice_idx)]) / dist281986) / sqrtl(dist281986)));
    }
    for (int j82036 = 0; j82036 < 4; j82036++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3708dinitial_mpg81946[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr3708dinitial_der82022[((24 * (0 + slice_idx)) + mapped_idx)] = -get_dfdx_cell(((j82036 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82023[(0 + slice_idx)] = (perturb_pos81911[((i82020 * 3) + slice_idx)] - (state81978[((j82036 * 6) + slice_idx)] + central_pos81909[(0 + slice_idx)]));
      }
      for (int loop_var82045 = 0; loop_var82045 < 24; loop_var82045++) {
        int mapped_idx;
        mapped_idx = loop_var82045;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81932[loop_var82045];
        ddist23094dinitial_der81985[mapped_idx] = ((((r82023[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022))) + ((r82023[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)))) + ((r82023[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022)) + (r82023[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022))));
      }
      dist281986 = (((r82023[0] * r82023[0]) + (r82023[1] * r82023[1])) + (r82023[2] * r82023[2]));
      for (int slice_idx = 0; slice_idx < (((j82036 * 3) + 3) - (j82036 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j82036 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81968[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81972[((24 * ((j82036 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j82036 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81969, dsat_acc2828dinitial_der81972) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81947, dr3708dinitial_der82022) * perturb_gms81910[i82020]) * dist281986) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985) * (r82023[(0 + slice_idx)] * perturb_gms81910[i82020]))) / (dist281986 * dist281986)) * sqrtl(dist281986)) - ((0.5 * (sqrtl((1.0 / dist281986)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985))) * ((r82023[(0 + slice_idx)] * perturb_gms81910[i82020]) / dist281986))) / (sqrtl(dist281986) * sqrtl(dist281986))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j82036 * 3) + 3) - (j82036 * 3)); slice_idx++) {
        sat_acc81973[((j82036 * 3) + slice_idx)] = (sat_acc81973[((j82036 * 3) + slice_idx)] + (((perturb_gms81910[i82020] * r82023[(0 + slice_idx)]) / dist281986) / sqrtl(dist281986)));
      }
    }
  }
  for (int i82051 = 1; i82051 < 4; i82051++) {
    long double r82054[3] = { 0.0 };
    long double dr4219dinitial_der82053[72] = { 0.0 };
    for (int j82055 = 0; j82055 < i82051; j82055++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr4219dinitial_mpg81952[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr4219dinitial_der82053[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j82055 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977) - get_dfdx_cell(((i82051 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82054[(0 + slice_idx)] = (state81978[((j82055 * 6) + slice_idx)] - state81978[((i82051 * 6) + slice_idx)]);
      }
      for (int loop_var82064 = 0; loop_var82064 < 24; loop_var82064++) {
        int mapped_idx;
        mapped_idx = loop_var82064;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81932[loop_var82064];
        ddist23094dinitial_der81985[mapped_idx] = ((((r82054[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053)) + (r82054[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053))) + ((r82054[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053)) + (r82054[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053)))) + ((r82054[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053)) + (r82054[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053))));
      }
      dist281986 = (((r82054[0] * r82054[0]) + (r82054[1] * r82054[1])) + (r82054[2] * r82054[2]));
      for (int slice_idx = 0; slice_idx < (((i82051 * 3) + 3) - (i82051 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((i82051 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81968[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81972[((24 * ((i82051 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82051 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81969, dsat_acc2828dinitial_der81972) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053) * sat_gms81912[j82055]) * dist281986) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985) * (r82054[(0 + slice_idx)] * sat_gms81912[j82055]))) / (dist281986 * dist281986)) * sqrtl(dist281986)) - ((0.5 * (sqrtl((1.0 / dist281986)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985))) * ((r82054[(0 + slice_idx)] * sat_gms81912[j82055]) / dist281986))) / (sqrtl(dist281986) * sqrtl(dist281986))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i82051 * 3) + 3) - (i82051 * 3)); slice_idx++) {
        sat_acc81973[((i82051 * 3) + slice_idx)] = (sat_acc81973[((i82051 * 3) + slice_idx)] + (((sat_gms81912[j82055] * r82054[(0 + slice_idx)]) / dist281986) / sqrtl(dist281986)));
      }
      for (int slice_idx = 0; slice_idx < (((j82055 * 3) + 3) - (j82055 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j82055 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81968[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81972[((24 * ((j82055 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j82055 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81969, dsat_acc2828dinitial_der81972) - (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81953, dr4219dinitial_der82053) * sat_gms81912[i82051]) * dist281986) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985) * (r82054[(0 + slice_idx)] * sat_gms81912[i82051]))) / (dist281986 * dist281986)) * sqrtl(dist281986)) - ((0.5 * (sqrtl((1.0 / dist281986)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81933, ddist23094dinitial_der81985))) * ((r82054[(0 + slice_idx)] * sat_gms81912[i82051]) / dist281986))) / (sqrtl(dist281986) * sqrtl(dist281986))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j82055 * 3) + 3) - (j82055 * 3)); slice_idx++) {
        sat_acc81973[((j82055 * 3) + slice_idx)] = (sat_acc81973[((j82055 * 3) + slice_idx)] - (((sat_gms81912[i82051] * r82054[(0 + slice_idx)]) / dist281986) / sqrtl(dist281986)));
      }
    }
  }
  long double grav_acc_local82075[3] = { 0.0 };
  long double dgrav_acc_local4544dinitial_der82074[72] = { 0.0 };
  long double grav_acc82077[3] = { 0.0 };
  long double dgrav_acc4569dinitial_der82076[72] = { 0.0 };
  long double acc82079[3] = { 0.0 };
  long double dacc4588dinitial_der82078[72] = { 0.0 };
  long double rot82080[9] = { 0.0 };
  long double v_082082[7] = { 0.0 };
  long double dv_04620dinitial_der82081[168] = { 0.0 };
  long double dv_0_dx82084[7] = { 0.0 };
  long double ddv_0_dx4638dinitial_der82083[168] = { 0.0 };
  long double dv_0_dy82086[7] = { 0.0 };
  long double ddv_0_dy4660dinitial_der82085[168] = { 0.0 };
  long double dv_0_dz82088[7] = { 0.0 };
  long double ddv_0_dz4682dinitial_der82087[168] = { 0.0 };
  {
    long double jupiter_rotation_matrix82097[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix82097, t81907);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot82080[(0 + slice_idx)] = jupiter_rotation_matrix82097[slice_idx];
    }
  }
  for (int i82098 = 0; i82098 < 4; i82098++) {
    long double x82101;
    long double dx4794dinitial_der82100[24] = { 0.0 };
    for (int loop_var82105 = 0; loop_var82105 < 24; loop_var82105++) {
      int mapped_idx;
      mapped_idx = loop_var82105;
      int al_index_name_symbol;
      al_index_name_symbol = dx4794dinitial_mpg81920[loop_var82105];
      dx4794dinitial_der82100[mapped_idx] = (get_dfdx_cell(((i82098 * 6) + 0), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977) / 0.0004778945025452157572e0);
    }
    x82101 = (state81978[((i82098 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y82108;
    long double dy4843dinitial_der82107[24] = { 0.0 };
    for (int loop_var82112 = 0; loop_var82112 < 24; loop_var82112++) {
      int mapped_idx;
      mapped_idx = loop_var82112;
      int al_index_name_symbol;
      al_index_name_symbol = dy4843dinitial_mpg81948[loop_var82112];
      dy4843dinitial_der82107[mapped_idx] = (get_dfdx_cell(((i82098 * 6) + 1), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977) / 0.0004778945025452157572e0);
    }
    y82108 = (state81978[((i82098 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z82115;
    long double dz4892dinitial_der82114[24] = { 0.0 };
    for (int loop_var82119 = 0; loop_var82119 < 24; loop_var82119++) {
      int mapped_idx;
      mapped_idx = loop_var82119;
      int al_index_name_symbol;
      al_index_name_symbol = dz4892dinitial_mpg81942[loop_var82119];
      dz4892dinitial_der82114[mapped_idx] = (get_dfdx_cell(((i82098 * 6) + 2), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81919, dstate2924dinitial_der81977) / 0.0004778945025452157572e0);
    }
    z82115 = (state81978[((i82098 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82122;
    long double d_x4942dinitial_der82121[24] = { 0.0 };
    for (int loop_var82126 = 0; loop_var82126 < 24; loop_var82126++) {
      int mapped_idx;
      mapped_idx = loop_var82126;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4942dinitial_mpg81954[loop_var82126];
      d_x4942dinitial_der82121[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81921, dx4794dinitial_der82100) * rot82080[0]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81949, dy4843dinitial_der82107) * rot82080[3])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81943, dz4892dinitial_der82114) * rot82080[6]));
    }
    _x82122 = (((rot82080[0] * x82101) + (rot82080[3] * y82108)) + (rot82080[6] * z82115));
    long double _y82129;
    long double d_y4993dinitial_der82128[24] = { 0.0 };
    for (int loop_var82133 = 0; loop_var82133 < 24; loop_var82133++) {
      int mapped_idx;
      mapped_idx = loop_var82133;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4993dinitial_mpg81922[loop_var82133];
      d_y4993dinitial_der82128[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81921, dx4794dinitial_der82100) * rot82080[1]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81949, dy4843dinitial_der82107) * rot82080[4])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81943, dz4892dinitial_der82114) * rot82080[7]));
    }
    _y82129 = (((rot82080[1] * x82101) + (rot82080[4] * y82108)) + (rot82080[7] * z82115));
    long double _z82136;
    long double d_z5044dinitial_der82135[24] = { 0.0 };
    for (int loop_var82140 = 0; loop_var82140 < 24; loop_var82140++) {
      int mapped_idx;
      mapped_idx = loop_var82140;
      int al_index_name_symbol;
      al_index_name_symbol = d_z5044dinitial_mpg81926[loop_var82140];
      d_z5044dinitial_der82135[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81921, dx4794dinitial_der82100) * rot82080[2]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81949, dy4843dinitial_der82107) * rot82080[5])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81943, dz4892dinitial_der82114) * rot82080[8]));
    }
    _z82136 = (((rot82080[2] * x82101) + (rot82080[5] * y82108)) + (rot82080[8] * z82115));
    long double _r282143;
    long double d_r25100dinitial_der82142[24] = { 0.0 };
    for (int loop_var82147 = 0; loop_var82147 < 24; loop_var82147++) {
      int mapped_idx;
      mapped_idx = loop_var82147;
      int al_index_name_symbol;
      al_index_name_symbol = d_r25100dinitial_mpg81938[loop_var82147];
      d_r25100dinitial_der82142[mapped_idx] = ((((_x82122 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121)) + (_x82122 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121))) + ((_y82129 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128)) + (_y82129 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128)))) + ((_z82136 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135))));
    }
    _r282143 = (((_x82122 * _x82122) + (_y82129 * _y82129)) + (_z82136 * _z82136));
    long double _r82150;
    long double d_r5143dinitial_der82149[24] = { 0.0 };
    for (int loop_var82154 = 0; loop_var82154 < 24; loop_var82154++) {
      int mapped_idx;
      mapped_idx = loop_var82154;
      int al_index_name_symbol;
      al_index_name_symbol = d_r5143dinitial_mpg81958[loop_var82154];
      d_r5143dinitial_der82149[mapped_idx] = (0.5 * (sqrtl((1.0 / _r282143)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142)));
    }
    _r82150 = sqrtl(_r282143);
    long double _r382157;
    long double d_r35167dinitial_der82156[24] = { 0.0 };
    for (int loop_var82161 = 0; loop_var82161 < 24; loop_var82161++) {
      int mapped_idx;
      mapped_idx = loop_var82161;
      int al_index_name_symbol;
      al_index_name_symbol = d_r35167dinitial_mpg81940[loop_var82161];
      d_r35167dinitial_der82156[mapped_idx] = ((_r82150 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142)) + (_r282143 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81959, d_r5143dinitial_der82149)));
    }
    _r382157 = (_r82150 * _r282143);
    long double _r482164;
    long double d_r45191dinitial_der82163[24] = { 0.0 };
    for (int loop_var82168 = 0; loop_var82168 < 24; loop_var82168++) {
      int mapped_idx;
      mapped_idx = loop_var82168;
      int al_index_name_symbol;
      al_index_name_symbol = d_r45191dinitial_mpg81966[loop_var82168];
      d_r45191dinitial_der82163[mapped_idx] = ((_r282143 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142)) + (_r282143 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142)));
    }
    _r482164 = (_r282143 * _r282143);
    long double _r582171;
    long double d_r55216dinitial_der82170[24] = { 0.0 };
    for (int loop_var82175 = 0; loop_var82175 < 24; loop_var82175++) {
      int mapped_idx;
      mapped_idx = loop_var82175;
      int al_index_name_symbol;
      al_index_name_symbol = d_r55216dinitial_mpg81914[loop_var82175];
      d_r55216dinitial_der82170[mapped_idx] = ((_r482164 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81959, d_r5143dinitial_der82149)) + (_r82150 * get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163)));
    }
    _r582171 = (_r482164 * _r82150);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81964[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der82081[((24 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r82150 * _r82150)) * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81959, d_r5143dinitial_der82149))));
        }
      }
    }
    v_082082[0] = (1.0e0 / _r82150);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81964[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der82081[((24 * 1) + mapped_idx)] = ((((get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * 1.7320508075688772936e0) * _r382157) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81941, d_r35167dinitial_der82156) * (_z82136 * 1.7320508075688772936e0))) / (_r382157 * _r382157));
        }
      }
    }
    v_082082[1] = ((1.7320508075688772936e0 * _z82136) / _r382157);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81928[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der82083[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121) * _r382157) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81941, d_r35167dinitial_der82156) * -_x82122)) / (_r382157 * _r382157));
        }
      }
    }
    dv_0_dx82084[0] = (-_x82122 / _r382157);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81960[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der82085[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128) * _r382157) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81941, d_r35167dinitial_der82156) * -_y82129)) / (_r382157 * _r382157));
        }
      }
    }
    dv_0_dy82086[0] = (-_y82129 / _r382157);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81934[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der82087[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * _r382157) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81941, d_r35167dinitial_der82156) * -_z82136)) / (_r382157 * _r382157));
        }
      }
    }
    dv_0_dz82088[0] = (-_z82136 / _r382157);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81928[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der82083[((24 * 1) + mapped_idx)] = ((((((_z82136 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121)) + (_x82122 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * -5.1961524227066318805e0))) * _r582171) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81915, d_r55216dinitial_der82170) * ((_z82136 * -5.1961524227066318805e0) * _x82122))) / (_r582171 * _r582171));
        }
      }
    }
    dv_0_dx82084[1] = (((-5.1961524227066318805e0 * _z82136) * _x82122) / _r582171);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81960[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der82085[((24 * 1) + mapped_idx)] = ((((((_z82136 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128)) + (_y82129 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * -5.1961524227066318805e0))) * _r582171) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81915, d_r55216dinitial_der82170) * ((_z82136 * -5.1961524227066318805e0) * _y82129))) / (_r582171 * _r582171));
        }
      }
    }
    dv_0_dy82086[1] = (((-5.1961524227066318805e0 * _z82136) * _y82129) / _r582171);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81934[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der82087[((24 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r382157 * _r382157)) * get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81941, d_r35167dinitial_der82156)))) - ((((((_z82136 * 3.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * 3.0e0))) * _r582171) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81915, d_r55216dinitial_der82170) * ((_z82136 * 3.0e0) * _z82136))) / (_r582171 * _r582171))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz82088[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382157) - (((3.0e0 * _z82136) * _z82136) / _r582171)));
    for (int n82209 = 2; n82209 < 7; n82209++) {
      long double coef182211;
      coef182211 = sqrtl(((((2.0e0 * ((long double) n82209)) - 1.0e0) * ((long double) ((2 * n82209) + 1))) / ((long double) (n82209 * n82209))));
      long double coef282215;
      coef282215 = sqrtl(((((((long double) n82209) - 1.0e0) * ((long double) (n82209 - 1))) * ((long double) ((2 * n82209) + 1))) / ((long double) ((n82209 * n82209) * ((2 * n82209) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82209));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_04620dinitial_mpg81964[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dv_04620dinitial_der82081[((24 * n82209) + mapped_idx)] = (((((((v_082082[(n82209 - 1)] * coef182211) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * (get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * coef182211))) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * ((v_082082[(n82209 - 1)] * coef182211) * _z82136))) / (_r282143 * _r282143)) - ((((get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * coef282215) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * (v_082082[(n82209 - 2)] * coef282215))) / (_r282143 * _r282143)));
          }
        }
      }
      v_082082[n82209] = ((((coef182211 * v_082082[(n82209 - 1)]) * _z82136) / _r282143) - ((coef282215 * v_082082[(n82209 - 2)]) / _r282143));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82209));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx4638dinitial_mpg81928[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dx4638dinitial_der82083[((24 * n82209) + mapped_idx)] = ((((_z82136 * coef182211) * ((((get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81929, ddv_0_dx4638dinitial_der82083) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * dv_0_dx82084[(n82209 - 1)])) / (_r282143 * _r282143)) - ((((((v_082082[(n82209 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121)) + (_x82122 * (get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((v_082082[(n82209 - 1)] * 2.0e0) * _x82122))) / (_r482164 * _r482164)))) + (((dv_0_dx82084[(n82209 - 1)] / _r282143) - (((v_082082[(n82209 - 1)] * 2.0e0) * _x82122) / _r482164)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * coef182211))) - (((((get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81929, ddv_0_dx4638dinitial_der82083) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * dv_0_dx82084[(n82209 - 2)])) / (_r282143 * _r282143)) - ((((((v_082082[(n82209 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81955, d_x4942dinitial_der82121)) + (_x82122 * (get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((v_082082[(n82209 - 2)] * 2.0e0) * _x82122))) / (_r482164 * _r482164))) * coef282215));
          }
        }
      }
      dv_0_dx82084[n82209] = (((coef182211 * _z82136) * ((dv_0_dx82084[(n82209 - 1)] / _r282143) - (((v_082082[(n82209 - 1)] * 2.0e0) * _x82122) / _r482164))) - (coef282215 * ((dv_0_dx82084[(n82209 - 2)] / _r282143) - (((v_082082[(n82209 - 2)] * 2.0e0) * _x82122) / _r482164))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82209));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy4660dinitial_mpg81960[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dy4660dinitial_der82085[((24 * n82209) + mapped_idx)] = ((((_z82136 * coef182211) * ((((get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81961, ddv_0_dy4660dinitial_der82085) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * dv_0_dy82086[(n82209 - 1)])) / (_r282143 * _r282143)) - ((((((v_082082[(n82209 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128)) + (_y82129 * (get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((v_082082[(n82209 - 1)] * 2.0e0) * _y82129))) / (_r482164 * _r482164)))) + (((dv_0_dy82086[(n82209 - 1)] / _r282143) - (((v_082082[(n82209 - 1)] * 2.0e0) * _y82129) / _r482164)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * coef182211))) - (((((get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81961, ddv_0_dy4660dinitial_der82085) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * dv_0_dy82086[(n82209 - 2)])) / (_r282143 * _r282143)) - ((((((v_082082[(n82209 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81923, d_y4993dinitial_der82128)) + (_y82129 * (get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((v_082082[(n82209 - 2)] * 2.0e0) * _y82129))) / (_r482164 * _r482164))) * coef282215));
          }
        }
      }
      dv_0_dy82086[n82209] = (((coef182211 * _z82136) * ((dv_0_dy82086[(n82209 - 1)] / _r282143) - (((v_082082[(n82209 - 1)] * 2.0e0) * _y82129) / _r482164))) - (coef282215 * ((dv_0_dy82086[(n82209 - 2)] / _r282143) - (((v_082082[(n82209 - 2)] * 2.0e0) * _y82129) / _r482164))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n82209));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz4682dinitial_mpg81934[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dz4682dinitial_der82087[((24 * n82209) + mapped_idx)] = ((((((((dv_0_dz82088[(n82209 - 1)] * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81935, ddv_0_dz4682dinitial_der82087))) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * (dv_0_dz82088[(n82209 - 1)] * _z82136))) / (_r282143 * _r282143)) + ((v_082082[(n82209 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r282143 * _r282143)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142)))) - ((((((_z82136 * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((_z82136 * 2.0e0) * _z82136))) / (_r482164 * _r482164)))) + (((1.0e0 / _r282143) - (((_z82136 * 2.0e0) * _z82136) / _r482164)) * get_dfdx_cell((n82209 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081)))) * coef182211) - (((((get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81935, ddv_0_dz4682dinitial_der82087) * _r282143) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81939, d_r25100dinitial_der82142) * dv_0_dz82088[(n82209 - 2)])) / (_r282143 * _r282143)) - ((((((v_082082[(n82209 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81927, d_z5044dinitial_der82135)) + (_z82136 * (get_dfdx_cell((n82209 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81965, dv_04620dinitial_der82081) * 2.0e0))) * _r482164) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81967, d_r45191dinitial_der82163) * ((v_082082[(n82209 - 2)] * 2.0e0) * _z82136))) / (_r482164 * _r482164))) * coef282215));
          }
        }
      }
      dv_0_dz82088[n82209] = ((coef182211 * (((dv_0_dz82088[(n82209 - 1)] * _z82136) / _r282143) + (v_082082[(n82209 - 1)] * ((1.0e0 / _r282143) - (((2.0e0 * _z82136) * _z82136) / _r482164))))) - (coef282215 * ((dv_0_dz82088[(n82209 - 2)] / _r282143) - (((v_082082[(n82209 - 2)] * 2.0e0) * _z82136) / _r482164))));
    }
    long double dpotential_dx82236;
    long double ddpotential_dx6387dinitial_der82235[24] = { 0.0 };
    for (int loop_var82240 = 0; loop_var82240 < 24; loop_var82240++) {
      int mapped_idx;
      mapped_idx = loop_var82240;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx6387dinitial_mpg81930[loop_var82240];
      ddpotential_dx6387dinitial_der82235[mapped_idx] = 0.0;
    }
    dpotential_dx82236 = ((long double) 0);
    long double dpotential_dy82242;
    long double ddpotential_dy6414dinitial_der82241[24] = { 0.0 };
    for (int loop_var82246 = 0; loop_var82246 < 24; loop_var82246++) {
      int mapped_idx;
      mapped_idx = loop_var82246;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy6414dinitial_mpg81970[loop_var82246];
      ddpotential_dy6414dinitial_der82241[mapped_idx] = 0.0;
    }
    dpotential_dy82242 = ((long double) 0);
    long double dpotential_dz82248;
    long double ddpotential_dz6441dinitial_der82247[24] = { 0.0 };
    for (int loop_var82252 = 0; loop_var82252 < 24; loop_var82252++) {
      int mapped_idx;
      mapped_idx = loop_var82252;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz6441dinitial_mpg81962[loop_var82252];
      ddpotential_dz6441dinitial_der82247[mapped_idx] = 0.0;
    }
    dpotential_dz82248 = ((long double) 0);
    for (int n82253 = 2; n82253 < 7; n82253++) {
      for (int loop_var82258 = 0; loop_var82258 < 24; loop_var82258++) {
        int mapped_idx;
        mapped_idx = loop_var82258;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx6387dinitial_mpg81930[loop_var82258];
        ddpotential_dx6387dinitial_der82235[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81931, ddpotential_dx6387dinitial_der82235) + (get_dfdx_cell(n82253, 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81929, ddv_0_dx4638dinitial_der82083) * central_grav81841[n82253]));
      }
      dpotential_dx82236 = (dpotential_dx82236 + (dv_0_dx82084[n82253] * central_grav81841[n82253]));
      for (int loop_var82263 = 0; loop_var82263 < 24; loop_var82263++) {
        int mapped_idx;
        mapped_idx = loop_var82263;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy6414dinitial_mpg81970[loop_var82263];
        ddpotential_dy6414dinitial_der82241[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81971, ddpotential_dy6414dinitial_der82241) + (get_dfdx_cell(n82253, 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81961, ddv_0_dy4660dinitial_der82085) * central_grav81841[n82253]));
      }
      dpotential_dy82242 = (dpotential_dy82242 + (dv_0_dy82086[n82253] * central_grav81841[n82253]));
      for (int loop_var82268 = 0; loop_var82268 < 24; loop_var82268++) {
        int mapped_idx;
        mapped_idx = loop_var82268;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz6441dinitial_mpg81962[loop_var82268];
        ddpotential_dz6441dinitial_der82247[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81963, ddpotential_dz6441dinitial_der82247) + (get_dfdx_cell(n82253, 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81935, ddv_0_dz4682dinitial_der82087) * central_grav81841[n82253]));
      }
      dpotential_dz82248 = (dpotential_dz82248 + (dv_0_dz82088[n82253] * central_grav81841[n82253]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81950[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der82074[((24 * 0) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81931, ddpotential_dx6387dinitial_der82235);
        }
      }
    }
    grav_acc_local82075[0] = dpotential_dx82236;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81950[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der82074[((24 * 1) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81971, ddpotential_dy6414dinitial_der82241);
        }
      }
    }
    grav_acc_local82075[1] = dpotential_dy82242;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 2));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81950[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der82074[((24 * 2) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81963, ddpotential_dz6441dinitial_der82247);
        }
      }
    }
    grav_acc_local82075[2] = dpotential_dz82248;
    for (int k82282 = 0; k82282 < 3; k82282++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * k82282));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc4569dinitial_mpg81956[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dgrav_acc4569dinitial_der82076[((24 * k82282) + mapped_idx)] = ((((get_dfdx_cell(0, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81951, dgrav_acc_local4544dinitial_der82074) * rot82080[((3 * k82282) + 0)]) + (get_dfdx_cell(1, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81951, dgrav_acc_local4544dinitial_der82074) * rot82080[((3 * k82282) + 1)])) + (get_dfdx_cell(2, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81951, dgrav_acc_local4544dinitial_der82074) * rot82080[((3 * k82282) + 2)])) / 2.2838315556293922983e-07);
          }
        }
      }
      grav_acc82077[k82282] = ((((rot82080[((3 * k82282) + 0)] * grav_acc_local82075[0]) + (rot82080[((3 * k82282) + 1)] * grav_acc_local82075[1])) + (rot82080[((3 * k82282) + 2)] * grav_acc_local82075[2])) / 2.2838315556293922983e-07);
    }
    for (int slice_idx = 0; slice_idx < (((i82098 * 3) + 3) - (i82098 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i82098 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81968[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81972[((24 * ((i82098 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82098 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81969, dsat_acc2828dinitial_der81972) + (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81957, dgrav_acc4569dinitial_der82076) * 2.8247609439046209905e-07));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i82098 * 3) + 3) - (i82098 * 3)); slice_idx++) {
      sat_acc81973[((i82098 * 3) + slice_idx)] = (sat_acc81973[((i82098 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc82077[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81936[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81974[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81937, dcentral_acc2846dinitial_der81974) - (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81957, dgrav_acc4569dinitial_der82076) * sat_gms81912[i82098]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81975[(0 + slice_idx)] = (central_acc81975[(0 + slice_idx)] - (sat_gms81912[i82098] * grav_acc82077[(0 + slice_idx)]));
    }
  }
  for (int i82296 = 0; i82296 < 4; i82296++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc4588dinitial_mpg81916[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dacc4588dinitial_der82078[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82296 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81969, dsat_acc2828dinitial_der81972) - get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81937, dcentral_acc2846dinitial_der81974));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc82079[(0 + slice_idx)] = (sat_acc81973[((i82296 * 3) + slice_idx)] - central_acc81975[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82296 * 6) + 3) - (i82296 * 6)); slice_idx++) {
        jupsatsystem81906[((i82296 * 6) + slice_idx)] = state81978[(((i82296 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82296 * 6) + 6) - ((i82296 * 6) + 3)); slice_idx++) {
        jupsatsystem81906[(((i82296 * 6) + 3) + slice_idx)] = acc82079[(0 + slice_idx)];
      }
    }
    for (int j82308 = 0; j82308 < 3; j82308++) {
      {
        for (int slice_idx = 0; slice_idx < ((24 + (((((i82296 * 6) + j82308) + 1) * 4) * 6)) - (24 + ((((i82296 * 6) + j82308) * 4) * 6))); slice_idx++) {
          jupsatsystem81906[((24 + ((((i82296 * 6) + j82308) * 4) * 6)) + slice_idx)] = state_and_derivatives81908[((24 + (((((i82296 * 6) + j82308) + 3) * 4) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((24 + (((((i82296 * 6) + j82308) + 4) * 4) * 6)) - (24 + (((((i82296 * 6) + j82308) + 3) * 4) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (24 + (((((i82296 * 6) + j82308) + 3) * 4) * 6)));
        jupsatsystem81906[func_slice_idx] = get_dfdx_cell(j82308, 24, al_index_name_symbol, 24, dacc4588dinitial_inv_mpg81917, dacc4588dinitial_der82078);
      }
    }
  }
  return 0;
}

int jupsatsystem_noderiv(long double *restrict jupsatsystem_noderiv82315, long double t82316, long double *restrict state82317, long double *restrict central_pos82318, long double *restrict perturb_gms82319, long double *restrict perturb_pos82320, long double *restrict sat_gms82321) {
  long double sat_acc82323[12] = { 0.0 };
  long double central_acc82324[3] = { 0.0 };
  long double dist282325;
  long double dist382326;
  for (int i82327 = 0; i82327 < 4; i82327++) {
    long double r82329[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82329[(0 + slice_idx)] = state82317[((i82327 * 6) + slice_idx)];
      }
    }
    dist282325 = (((r82329[0] * r82329[0]) + (r82329[1] * r82329[1])) + (r82329[2] * r82329[2]));
    dist382326 = (dist282325 * sqrtl(dist282325));
    {
      for (int slice_idx = 0; slice_idx < (((i82327 * 3) + 3) - (i82327 * 3)); slice_idx++) {
        sat_acc82323[((i82327 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r82329[(0 + slice_idx)]) / dist382326);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82324[(0 + slice_idx)] = (central_acc82324[(0 + slice_idx)] + ((sat_gms82321[i82327] * r82329[(0 + slice_idx)]) / dist382326));
      }
    }
  }
  for (int i82345 = 0; i82345 < 4; i82345++) {
    long double r82347[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82347[(0 + slice_idx)] = (perturb_pos82320[((i82345 * 3) + slice_idx)] - central_pos82318[(0 + slice_idx)]);
      }
    }
    dist282325 = (((r82347[0] * r82347[0]) + (r82347[1] * r82347[1])) + (r82347[2] * r82347[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82324[(0 + slice_idx)] = (central_acc82324[(0 + slice_idx)] + (((perturb_gms82319[i82345] * r82347[(0 + slice_idx)]) / dist282325) / sqrtl(dist282325)));
      }
    }
    for (int j82357 = 0; j82357 < 4; j82357++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82347[(0 + slice_idx)] = (perturb_pos82320[((i82345 * 3) + slice_idx)] - (state82317[((j82357 * 6) + slice_idx)] + central_pos82318[(0 + slice_idx)]));
        }
      }
      dist282325 = (((r82347[0] * r82347[0]) + (r82347[1] * r82347[1])) + (r82347[2] * r82347[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j82357 * 3) + 3) - (j82357 * 3)); slice_idx++) {
          sat_acc82323[((j82357 * 3) + slice_idx)] = (sat_acc82323[((j82357 * 3) + slice_idx)] + (((perturb_gms82319[i82345] * r82347[(0 + slice_idx)]) / dist282325) / sqrtl(dist282325)));
        }
      }
    }
  }
  for (int i82368 = 1; i82368 < 4; i82368++) {
    long double r82370[3] = { 0.0 };
    for (int j82371 = 0; j82371 < i82368; j82371++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82370[(0 + slice_idx)] = (state82317[((j82371 * 6) + slice_idx)] - state82317[((i82368 * 6) + slice_idx)]);
        }
      }
      dist282325 = (((r82370[0] * r82370[0]) + (r82370[1] * r82370[1])) + (r82370[2] * r82370[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i82368 * 3) + 3) - (i82368 * 3)); slice_idx++) {
          sat_acc82323[((i82368 * 3) + slice_idx)] = (sat_acc82323[((i82368 * 3) + slice_idx)] + (((sat_gms82321[j82371] * r82370[(0 + slice_idx)]) / dist282325) / sqrtl(dist282325)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j82371 * 3) + 3) - (j82371 * 3)); slice_idx++) {
          sat_acc82323[((j82371 * 3) + slice_idx)] = (sat_acc82323[((j82371 * 3) + slice_idx)] - (((sat_gms82321[i82368] * r82370[(0 + slice_idx)]) / dist282325) / sqrtl(dist282325)));
        }
      }
    }
  }
  long double grav_acc_local82385[3] = { 0.0 };
  long double grav_acc82386[3] = { 0.0 };
  long double acc82387[3] = { 0.0 };
  long double rot82388[9] = { 0.0 };
  long double v_082389[7] = { 0.0 };
  long double dv_0_dx82390[7] = { 0.0 };
  long double dv_0_dy82391[7] = { 0.0 };
  long double dv_0_dz82392[7] = { 0.0 };
  {
    long double jupiter_rotation_matrix82401[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix82401, t82316);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot82388[(0 + slice_idx)] = jupiter_rotation_matrix82401[slice_idx];
    }
  }
  for (int i82402 = 0; i82402 < 4; i82402++) {
    long double x82404;
    x82404 = (state82317[((i82402 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y82408;
    y82408 = (state82317[((i82402 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z82412;
    z82412 = (state82317[((i82402 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82416;
    _x82416 = (((rot82388[0] * x82404) + (rot82388[3] * y82408)) + (rot82388[6] * z82412));
    long double _y82420;
    _y82420 = (((rot82388[1] * x82404) + (rot82388[4] * y82408)) + (rot82388[7] * z82412));
    long double _z82424;
    _z82424 = (((rot82388[2] * x82404) + (rot82388[5] * y82408)) + (rot82388[8] * z82412));
    long double _r282428;
    _r282428 = (((_x82416 * _x82416) + (_y82420 * _y82420)) + (_z82424 * _z82424));
    long double r82432;
    r82432 = sqrtl(_r282428);
    long double _r382436;
    _r382436 = (r82432 * _r282428);
    long double _r482440;
    _r482440 = (_r282428 * _r282428);
    long double _r582444;
    _r582444 = (_r482440 * r82432);
    v_082389[0] = (1.0e0 / r82432);
    v_082389[1] = ((1.7320508075688772936e0 * _z82424) / _r382436);
    dv_0_dx82390[0] = (-_x82416 / _r382436);
    dv_0_dy82391[0] = (-_y82420 / _r382436);
    dv_0_dz82392[0] = (-_z82424 / _r382436);
    dv_0_dx82390[1] = (((-5.1961524227066318805e0 * _z82424) * _x82416) / _r582444);
    dv_0_dy82391[1] = (((-5.1961524227066318805e0 * _z82424) * _y82420) / _r582444);
    dv_0_dz82392[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382436) - (((3.0e0 * _z82424) * _z82424) / _r582444)));
    for (int n82472 = 2; n82472 < 7; n82472++) {
      long double coef182474;
      coef182474 = sqrtl(((((2.0e0 * ((long double) n82472)) - 1.0e0) * ((long double) ((2 * n82472) + 1))) / ((long double) (n82472 * n82472))));
      long double coef282478;
      coef282478 = sqrtl(((((((long double) n82472) - 1.0e0) * ((long double) (n82472 - 1))) * ((long double) ((2 * n82472) + 1))) / ((long double) ((n82472 * n82472) * ((2 * n82472) - 3)))));
      v_082389[n82472] = ((((coef182474 * v_082389[(n82472 - 1)]) * _z82424) / _r282428) - ((coef282478 * v_082389[(n82472 - 2)]) / _r282428));
      dv_0_dx82390[n82472] = (((coef182474 * _z82424) * ((dv_0_dx82390[(n82472 - 1)] / _r282428) - (((v_082389[(n82472 - 1)] * 2.0e0) * _x82416) / _r482440))) - (coef282478 * ((dv_0_dx82390[(n82472 - 2)] / _r282428) - (((v_082389[(n82472 - 2)] * 2.0e0) * _x82416) / _r482440))));
      dv_0_dy82391[n82472] = (((coef182474 * _z82424) * ((dv_0_dy82391[(n82472 - 1)] / _r282428) - (((v_082389[(n82472 - 1)] * 2.0e0) * _y82420) / _r482440))) - (coef282478 * ((dv_0_dy82391[(n82472 - 2)] / _r282428) - (((v_082389[(n82472 - 2)] * 2.0e0) * _y82420) / _r482440))));
      dv_0_dz82392[n82472] = ((coef182474 * (((dv_0_dz82392[(n82472 - 1)] * _z82424) / _r282428) + (v_082389[(n82472 - 1)] * ((1.0e0 / _r282428) - (((2.0e0 * _z82424) * _z82424) / _r482440))))) - (coef282478 * ((dv_0_dz82392[(n82472 - 2)] / _r282428) - (((v_082389[(n82472 - 2)] * 2.0e0) * _z82424) / _r482440))));
    }
    long double dpotential_dx82494;
    dpotential_dx82494 = ((long double) 0);
    long double dpotential_dy82498;
    dpotential_dy82498 = ((long double) 0);
    long double dpotential_dz82502;
    dpotential_dz82502 = ((long double) 0);
    for (int n82506 = 2; n82506 < 7; n82506++) {
      dpotential_dx82494 = (dpotential_dx82494 + (dv_0_dx82390[n82506] * central_grav81841[n82506]));
      dpotential_dy82498 = (dpotential_dy82498 + (dv_0_dy82391[n82506] * central_grav81841[n82506]));
      dpotential_dz82502 = (dpotential_dz82502 + (dv_0_dz82392[n82506] * central_grav81841[n82506]));
    }
    grav_acc_local82385[0] = dpotential_dx82494;
    grav_acc_local82385[1] = dpotential_dy82498;
    grav_acc_local82385[2] = dpotential_dz82502;
    for (int k82526 = 0; k82526 < 3; k82526++) {
      grav_acc82386[k82526] = ((((rot82388[((3 * k82526) + 0)] * grav_acc_local82385[0]) + (rot82388[((3 * k82526) + 1)] * grav_acc_local82385[1])) + (rot82388[((3 * k82526) + 2)] * grav_acc_local82385[2])) / 2.2838315556293922983e-07);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82402 * 3) + 3) - (i82402 * 3)); slice_idx++) {
        sat_acc82323[((i82402 * 3) + slice_idx)] = (sat_acc82323[((i82402 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc82386[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82324[(0 + slice_idx)] = (central_acc82324[(0 + slice_idx)] - (sat_gms82321[i82402] * grav_acc82386[(0 + slice_idx)]));
      }
    }
  }
  for (int i82537 = 0; i82537 < 4; i82537++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc82387[(0 + slice_idx)] = (sat_acc82323[((i82537 * 3) + slice_idx)] - central_acc82324[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82537 * 6) + 3) - (i82537 * 6)); slice_idx++) {
        jupsatsystem_noderiv82315[((i82537 * 6) + slice_idx)] = state82317[(((i82537 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82537 * 6) + 6) - ((i82537 * 6) + 3)); slice_idx++) {
        jupsatsystem_noderiv82315[(((i82537 * 6) + 3) + slice_idx)] = acc82387[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

