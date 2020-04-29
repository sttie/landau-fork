#include <math.h>

static inline long double get_dfdx_cell(int full_idx81555, int dx_mapped_size81556, int al_idx81557, int inv_mapping_period81558, int *restrict dx_idx_mappings81559, long double *restrict der_vec81560) {
  return ((al_idx81557 < inv_mapping_period81558) ? ((dx_idx_mappings81559[((inv_mapping_period81558 * full_idx81555) + al_idx81557)] >= 0) ? der_vec81560[((dx_mapped_size81556 * full_idx81555) + dx_idx_mappings81559[((inv_mapping_period81558 * full_idx81555) + al_idx81557)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx81555, int dx_mapped_size81556, int al_idx81557, int inv_mapping_period81558, int *restrict dx_idx_mappings81559, long double *restrict der_vec81560) {
  return ((al_idx81557 < inv_mapping_period81558) ? ((al_idx81557 == full_idx81555) ? 1.0 : ((dx_idx_mappings81559[((inv_mapping_period81558 * full_idx81555) + al_idx81557)] >= 0) ? der_vec81560[((dx_mapped_size81556 * full_idx81555) + dx_idx_mappings81559[((inv_mapping_period81558 * full_idx81555) + al_idx81557)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx81562, int inv_mapping_period81563, int *restrict dx_idx_mappings81564, long double *restrict der_vec81565) {
  return ((((al_idx81562 < inv_mapping_period81563) ? dx_idx_mappings81564[al_idx81562] : -1) >= 0) ? der_vec81565[((al_idx81562 < inv_mapping_period81563) ? dx_idx_mappings81564[al_idx81562] : -1)] : 0.0);
}

const long double pole_ra81566[2] = { 268.05659500000001572e0, -0.006498999999999999569e0 };
const long double pole_dec81567[2] = { 64.49530300000000693e0, 0.0024130000000000002142e0 };
const long double pm81568[2] = { 284.94999999999998863e0, 870.5359999999999445e0 };
const long double nut_prec_ra81569[5] = { 0.00011699999999999999788e0, 0.00093800000000000003167e0, 0.0014319999999999998962e0, 3.000000000000000076e-05, 0.0021500000000000000014e0 };
const long double nut_prec_dec81570[5] = { 5.0000000000000002396e-05, 0.00040400000000000000798e0, 0.0006170000000000000354e0, -1.29999999999999992e-05, 0.0009259999999999999568e0 };
const long double Jabcde_081571[5] = { 99.36071400000000153e0, 175.89536899999998809e0, 300.32316200000002482e0, 114.01230499999999779e0, 49.511251000000001454e0 };
const long double Jabcde_T81572[5] = { 4850.4045999999998457e0, 1191.9604999999999109e0, 262.54750000000001364e0, 6070.247599999999693e0, 64.29999999999999716e0 };
const long double central_grav81573[7] = { 0.0e0, 0.0e0, -0.0065724808672554692115e0, 0.0e0, 0.00019554099999999999462e0, 0.0e0, -9.4975767597683728686e-06 };

int jupiter_rotation_matrix(long double *restrict jupiter_rotation_matrix81574, long double t81575) {
  long double T81577;
  T81577 = (t81575 / 36525.0e0);
  long double alpha_081581;
  alpha_081581 = (268.05659500000001572e0 + (-0.006498999999999999569e0 * T81577));
  long double delta_081585;
  delta_081585 = (64.49530300000000693e0 + (0.0024130000000000002142e0 * T81577));
  long double W81589;
  W81589 = (((284.94999999999998863e0 + (870.5359999999999445e0 * t81575)) * 3.141592653589793116e0) / 180.0e0);
  for (int i81593 = 0; i81593 < 5; i81593++) {
    long double J81595;
    J81595 = (((Jabcde_081571[i81593] + (Jabcde_T81572[i81593] * T81577)) * 3.141592653589793116e0) / 180.0e0);
    alpha_081581 = (alpha_081581 + (nut_prec_ra81569[i81593] * sinl(J81595)));
    delta_081585 = (delta_081585 + (nut_prec_dec81570[i81593] * cosl(J81595)));
  }
  alpha_081581 = (alpha_081581 * 0.017453292519943295088e0);
  delta_081585 = (delta_081585 * 0.017453292519943295088e0);
  jupiter_rotation_matrix81574[0] = ((-sinl(alpha_081581) * cosl(W81589)) - ((cosl(alpha_081581) * sinl(delta_081585)) * sinl(W81589)));
  jupiter_rotation_matrix81574[1] = ((sinl(alpha_081581) * sinl(W81589)) - ((cosl(alpha_081581) * sinl(delta_081585)) * cosl(W81589)));
  jupiter_rotation_matrix81574[2] = (cosl(alpha_081581) * cosl(delta_081585));
  jupiter_rotation_matrix81574[3] = ((cosl(alpha_081581) * cosl(W81589)) - ((sinl(alpha_081581) * sinl(delta_081585)) * sinl(W81589)));
  jupiter_rotation_matrix81574[4] = ((-cosl(alpha_081581) * sinl(W81589)) - ((sinl(alpha_081581) * sinl(delta_081585)) * cosl(W81589)));
  jupiter_rotation_matrix81574[5] = (cosl(delta_081585) * sinl(alpha_081581));
  jupiter_rotation_matrix81574[6] = (cosl(delta_081585) * sinl(W81589));
  jupiter_rotation_matrix81574[7] = (cosl(delta_081585) * cosl(W81589));
  jupiter_rotation_matrix81574[8] = sinl(delta_081585);
  return 0;
}

int jupsatsystem(long double *restrict jupsatsystem81638, long double t81639, long double *restrict state_and_derivatives81640, long double *restrict central_pos81641, long double *restrict perturb_gms81642, long double *restrict perturb_pos81643, long double *restrict sat_gms81644) {
  static int d_r55216dinitial_mpg81646[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_inv_mpg81647[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_mpg81648[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_inv_mpg81649[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dstate2924dinitial_mpg81650[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2924dinitial_inv_mpg81651[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dx4794dinitial_mpg81652[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_inv_mpg81653[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_mpg81654[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_inv_mpg81655[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_mpg81656[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_inv_mpg81657[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_mpg81658[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_inv_mpg81659[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_mpg81660[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_inv_mpg81661[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_mpg81662[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_inv_mpg81663[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_mpg81664[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_inv_mpg81665[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_mpg81666[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_inv_mpg81667[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_mpg81668[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_inv_mpg81669[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_mpg81670[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_inv_mpg81671[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_mpg81672[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_inv_mpg81673[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_mpg81674[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_inv_mpg81675[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_mpg81676[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_inv_mpg81677[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_mpg81678[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_inv_mpg81679[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_mpg81680[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_inv_mpg81681[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_mpg81682[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_inv_mpg81683[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_mpg81684[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_inv_mpg81685[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_mpg81686[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_inv_mpg81687[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_mpg81688[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_inv_mpg81689[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_mpg81690[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_inv_mpg81691[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_mpg81692[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_inv_mpg81693[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_mpg81694[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_inv_mpg81695[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dv_04620dinitial_mpg81696[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_04620dinitial_inv_mpg81697[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int d_r45191dinitial_mpg81698[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_inv_mpg81699[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_mpg81700[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_inv_mpg81701[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_mpg81702[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_inv_mpg81703[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  long double sat_acc81705[12] = { 0.0 };
  long double dsat_acc2828dinitial_der81704[288] = { 0.0 };
  long double central_acc81707[3] = { 0.0 };
  long double dcentral_acc2846dinitial_der81706[72] = { 0.0 };
  long double state_derivatives_initial81708[576] = { 0.0 };
  long double state81710[24] = { 0.0 };
  long double dstate2924dinitial_der81709[576] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 576; slice_idx++) {
      state_derivatives_initial81708[(0 + slice_idx)] = state_and_derivatives81640[(24 + slice_idx)];
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
        al_index_name_symbol = dstate2924dinitial_mpg81650[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dstate2924dinitial_der81709[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    state81710[(0 + slice_idx)] = state_and_derivatives81640[(0 + slice_idx)];
  }
  long double dist281718;
  long double ddist23094dinitial_der81717[24] = { 0.0 };
  long double dist381720;
  long double ddist33107dinitial_der81719[24] = { 0.0 };
  for (int i81721 = 0; i81721 < 24; i81721++) {
    for (int j81723 = 0; j81723 < 24; j81723++) {
      {
        int mapped_idx;
        mapped_idx = ((j81723 < 24) ? dstate2924dinitial_inv_mpg81651[((i81721 * 24) + j81723)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2924dinitial_der81709[((i81721 * 24) + mapped_idx)] = state_derivatives_initial81708[(((i81721 * 4) * 6) + j81723)];
        }
      }
    }
  }
  for (int i81726 = 0; i81726 < 4; i81726++) {
    long double r81729[3] = { 0.0 };
    long double dr3376dinitial_der81728[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3376dinitial_mpg81676[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3376dinitial_der81728[((24 * (0 + slice_idx)) + mapped_idx)] = get_dfdx_cell(((i81726 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81729[(0 + slice_idx)] = state81710[((i81726 * 6) + slice_idx)];
    }
    for (int loop_var81737 = 0; loop_var81737 < 24; loop_var81737++) {
      int mapped_idx;
      mapped_idx = loop_var81737;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81664[loop_var81737];
      ddist23094dinitial_der81717[mapped_idx] = ((((r81729[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728)) + (r81729[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728))) + ((r81729[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728)) + (r81729[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728)))) + ((r81729[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728)) + (r81729[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728))));
    }
    dist281718 = (((r81729[0] * r81729[0]) + (r81729[1] * r81729[1])) + (r81729[2] * r81729[2]));
    for (int loop_var81742 = 0; loop_var81742 < 24; loop_var81742++) {
      int mapped_idx;
      mapped_idx = loop_var81742;
      int al_index_name_symbol;
      al_index_name_symbol = ddist33107dinitial_mpg81656[loop_var81742];
      ddist33107dinitial_der81719[mapped_idx] = ((dist281718 * (0.5 * (sqrtl((1.0 / dist281718)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717)))) + (sqrtl(dist281718) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717)));
    }
    dist381720 = (dist281718 * sqrtl(dist281718));
    for (int slice_idx = 0; slice_idx < (((i81726 * 3) + 3) - (i81726 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81726 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81700[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81704[((24 * ((i81726 * 3) + slice_idx)) + mapped_idx)] = ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728) * -2.8247609439046209905e-07) * dist381720) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81657, ddist33107dinitial_der81719) * (r81729[(0 + slice_idx)] * -2.8247609439046209905e-07))) / (dist381720 * dist381720));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81726 * 3) + 3) - (i81726 * 3)); slice_idx++) {
      sat_acc81705[((i81726 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r81729[(0 + slice_idx)]) / dist381720);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81668[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81706[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81669, dcentral_acc2846dinitial_der81706) + ((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3376dinitial_inv_mpg81677, dr3376dinitial_der81728) * sat_gms81644[i81726]) * dist381720) - (get_dfdx_var(al_index_name_symbol, 24, ddist33107dinitial_inv_mpg81657, ddist33107dinitial_der81719) * (r81729[(0 + slice_idx)] * sat_gms81644[i81726]))) / (dist381720 * dist381720)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81707[(0 + slice_idx)] = (central_acc81707[(0 + slice_idx)] + ((sat_gms81644[i81726] * r81729[(0 + slice_idx)]) / dist381720));
    }
  }
  for (int i81752 = 0; i81752 < 4; i81752++) {
    long double r81755[3] = { 0.0 };
    long double dr3708dinitial_der81754[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3708dinitial_mpg81678[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3708dinitial_der81754[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r81755[(0 + slice_idx)] = (perturb_pos81643[((i81752 * 3) + slice_idx)] - central_pos81641[(0 + slice_idx)]);
    }
    for (int loop_var81762 = 0; loop_var81762 < 24; loop_var81762++) {
      int mapped_idx;
      mapped_idx = loop_var81762;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg81664[loop_var81762];
      ddist23094dinitial_der81717[mapped_idx] = ((((r81755[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754))) + ((r81755[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)))) + ((r81755[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754))));
    }
    dist281718 = (((r81755[0] * r81755[0]) + (r81755[1] * r81755[1])) + (r81755[2] * r81755[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81668[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81706[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81669, dcentral_acc2846dinitial_der81706) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754) * perturb_gms81642[i81752]) * dist281718) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717) * (r81755[(0 + slice_idx)] * perturb_gms81642[i81752]))) / (dist281718 * dist281718)) * sqrtl(dist281718)) - ((0.5 * (sqrtl((1.0 / dist281718)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717))) * ((r81755[(0 + slice_idx)] * perturb_gms81642[i81752]) / dist281718))) / (sqrtl(dist281718) * sqrtl(dist281718))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81707[(0 + slice_idx)] = (central_acc81707[(0 + slice_idx)] + (((perturb_gms81642[i81752] * r81755[(0 + slice_idx)]) / dist281718) / sqrtl(dist281718)));
    }
    for (int j81768 = 0; j81768 < 4; j81768++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3708dinitial_mpg81678[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr3708dinitial_der81754[((24 * (0 + slice_idx)) + mapped_idx)] = -get_dfdx_cell(((j81768 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81755[(0 + slice_idx)] = (perturb_pos81643[((i81752 * 3) + slice_idx)] - (state81710[((j81768 * 6) + slice_idx)] + central_pos81641[(0 + slice_idx)]));
      }
      for (int loop_var81777 = 0; loop_var81777 < 24; loop_var81777++) {
        int mapped_idx;
        mapped_idx = loop_var81777;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81664[loop_var81777];
        ddist23094dinitial_der81717[mapped_idx] = ((((r81755[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754))) + ((r81755[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)))) + ((r81755[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754)) + (r81755[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754))));
      }
      dist281718 = (((r81755[0] * r81755[0]) + (r81755[1] * r81755[1])) + (r81755[2] * r81755[2]));
      for (int slice_idx = 0; slice_idx < (((j81768 * 3) + 3) - (j81768 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81768 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81700[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81704[((24 * ((j81768 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81768 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81701, dsat_acc2828dinitial_der81704) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr3708dinitial_inv_mpg81679, dr3708dinitial_der81754) * perturb_gms81642[i81752]) * dist281718) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717) * (r81755[(0 + slice_idx)] * perturb_gms81642[i81752]))) / (dist281718 * dist281718)) * sqrtl(dist281718)) - ((0.5 * (sqrtl((1.0 / dist281718)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717))) * ((r81755[(0 + slice_idx)] * perturb_gms81642[i81752]) / dist281718))) / (sqrtl(dist281718) * sqrtl(dist281718))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81768 * 3) + 3) - (j81768 * 3)); slice_idx++) {
        sat_acc81705[((j81768 * 3) + slice_idx)] = (sat_acc81705[((j81768 * 3) + slice_idx)] + (((perturb_gms81642[i81752] * r81755[(0 + slice_idx)]) / dist281718) / sqrtl(dist281718)));
      }
    }
  }
  for (int i81783 = 1; i81783 < 4; i81783++) {
    long double r81786[3] = { 0.0 };
    long double dr4219dinitial_der81785[72] = { 0.0 };
    for (int j81787 = 0; j81787 < i81783; j81787++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr4219dinitial_mpg81684[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr4219dinitial_der81785[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81787 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709) - get_dfdx_cell(((i81783 * 6) + slice_idx), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r81786[(0 + slice_idx)] = (state81710[((j81787 * 6) + slice_idx)] - state81710[((i81783 * 6) + slice_idx)]);
      }
      for (int loop_var81796 = 0; loop_var81796 < 24; loop_var81796++) {
        int mapped_idx;
        mapped_idx = loop_var81796;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg81664[loop_var81796];
        ddist23094dinitial_der81717[mapped_idx] = ((((r81786[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785)) + (r81786[0] * get_dfdx_cell(0, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785))) + ((r81786[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785)) + (r81786[1] * get_dfdx_cell(1, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785)))) + ((r81786[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785)) + (r81786[2] * get_dfdx_cell(2, 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785))));
      }
      dist281718 = (((r81786[0] * r81786[0]) + (r81786[1] * r81786[1])) + (r81786[2] * r81786[2]));
      for (int slice_idx = 0; slice_idx < (((i81783 * 3) + 3) - (i81783 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((i81783 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81700[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81704[((24 * ((i81783 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81783 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81701, dsat_acc2828dinitial_der81704) + (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785) * sat_gms81644[j81787]) * dist281718) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717) * (r81786[(0 + slice_idx)] * sat_gms81644[j81787]))) / (dist281718 * dist281718)) * sqrtl(dist281718)) - ((0.5 * (sqrtl((1.0 / dist281718)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717))) * ((r81786[(0 + slice_idx)] * sat_gms81644[j81787]) / dist281718))) / (sqrtl(dist281718) * sqrtl(dist281718))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i81783 * 3) + 3) - (i81783 * 3)); slice_idx++) {
        sat_acc81705[((i81783 * 3) + slice_idx)] = (sat_acc81705[((i81783 * 3) + slice_idx)] + (((sat_gms81644[j81787] * r81786[(0 + slice_idx)]) / dist281718) / sqrtl(dist281718)));
      }
      for (int slice_idx = 0; slice_idx < (((j81787 * 3) + 3) - (j81787 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j81787 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg81700[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der81704[((24 * ((j81787 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((j81787 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81701, dsat_acc2828dinitial_der81704) - (((((((get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dr4219dinitial_inv_mpg81685, dr4219dinitial_der81785) * sat_gms81644[i81783]) * dist281718) - (get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717) * (r81786[(0 + slice_idx)] * sat_gms81644[i81783]))) / (dist281718 * dist281718)) * sqrtl(dist281718)) - ((0.5 * (sqrtl((1.0 / dist281718)) * get_dfdx_var(al_index_name_symbol, 24, ddist23094dinitial_inv_mpg81665, ddist23094dinitial_der81717))) * ((r81786[(0 + slice_idx)] * sat_gms81644[i81783]) / dist281718))) / (sqrtl(dist281718) * sqrtl(dist281718))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j81787 * 3) + 3) - (j81787 * 3)); slice_idx++) {
        sat_acc81705[((j81787 * 3) + slice_idx)] = (sat_acc81705[((j81787 * 3) + slice_idx)] - (((sat_gms81644[i81783] * r81786[(0 + slice_idx)]) / dist281718) / sqrtl(dist281718)));
      }
    }
  }
  long double grav_acc_local81807[3] = { 0.0 };
  long double dgrav_acc_local4544dinitial_der81806[72] = { 0.0 };
  long double grav_acc81809[3] = { 0.0 };
  long double dgrav_acc4569dinitial_der81808[72] = { 0.0 };
  long double acc81811[3] = { 0.0 };
  long double dacc4588dinitial_der81810[72] = { 0.0 };
  long double rot81812[9] = { 0.0 };
  long double v_081814[7] = { 0.0 };
  long double dv_04620dinitial_der81813[168] = { 0.0 };
  long double dv_0_dx81816[7] = { 0.0 };
  long double ddv_0_dx4638dinitial_der81815[168] = { 0.0 };
  long double dv_0_dy81818[7] = { 0.0 };
  long double ddv_0_dy4660dinitial_der81817[168] = { 0.0 };
  long double dv_0_dz81820[7] = { 0.0 };
  long double ddv_0_dz4682dinitial_der81819[168] = { 0.0 };
  {
    long double jupiter_rotation_matrix81829[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix81829, t81639);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot81812[(0 + slice_idx)] = jupiter_rotation_matrix81829[slice_idx];
    }
  }
  for (int i81830 = 0; i81830 < 4; i81830++) {
    long double x81833;
    long double dx4794dinitial_der81832[24] = { 0.0 };
    for (int loop_var81837 = 0; loop_var81837 < 24; loop_var81837++) {
      int mapped_idx;
      mapped_idx = loop_var81837;
      int al_index_name_symbol;
      al_index_name_symbol = dx4794dinitial_mpg81652[loop_var81837];
      dx4794dinitial_der81832[mapped_idx] = (get_dfdx_cell(((i81830 * 6) + 0), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709) / 0.0004778945025452157572e0);
    }
    x81833 = (state81710[((i81830 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y81840;
    long double dy4843dinitial_der81839[24] = { 0.0 };
    for (int loop_var81844 = 0; loop_var81844 < 24; loop_var81844++) {
      int mapped_idx;
      mapped_idx = loop_var81844;
      int al_index_name_symbol;
      al_index_name_symbol = dy4843dinitial_mpg81680[loop_var81844];
      dy4843dinitial_der81839[mapped_idx] = (get_dfdx_cell(((i81830 * 6) + 1), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709) / 0.0004778945025452157572e0);
    }
    y81840 = (state81710[((i81830 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z81847;
    long double dz4892dinitial_der81846[24] = { 0.0 };
    for (int loop_var81851 = 0; loop_var81851 < 24; loop_var81851++) {
      int mapped_idx;
      mapped_idx = loop_var81851;
      int al_index_name_symbol;
      al_index_name_symbol = dz4892dinitial_mpg81674[loop_var81851];
      dz4892dinitial_der81846[mapped_idx] = (get_dfdx_cell(((i81830 * 6) + 2), 24, al_index_name_symbol, 24, dstate2924dinitial_inv_mpg81651, dstate2924dinitial_der81709) / 0.0004778945025452157572e0);
    }
    z81847 = (state81710[((i81830 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x81854;
    long double d_x4942dinitial_der81853[24] = { 0.0 };
    for (int loop_var81858 = 0; loop_var81858 < 24; loop_var81858++) {
      int mapped_idx;
      mapped_idx = loop_var81858;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4942dinitial_mpg81686[loop_var81858];
      d_x4942dinitial_der81853[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81653, dx4794dinitial_der81832) * rot81812[0]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81681, dy4843dinitial_der81839) * rot81812[3])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81675, dz4892dinitial_der81846) * rot81812[6]));
    }
    _x81854 = (((rot81812[0] * x81833) + (rot81812[3] * y81840)) + (rot81812[6] * z81847));
    long double _y81861;
    long double d_y4993dinitial_der81860[24] = { 0.0 };
    for (int loop_var81865 = 0; loop_var81865 < 24; loop_var81865++) {
      int mapped_idx;
      mapped_idx = loop_var81865;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4993dinitial_mpg81654[loop_var81865];
      d_y4993dinitial_der81860[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81653, dx4794dinitial_der81832) * rot81812[1]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81681, dy4843dinitial_der81839) * rot81812[4])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81675, dz4892dinitial_der81846) * rot81812[7]));
    }
    _y81861 = (((rot81812[1] * x81833) + (rot81812[4] * y81840)) + (rot81812[7] * z81847));
    long double _z81868;
    long double d_z5044dinitial_der81867[24] = { 0.0 };
    for (int loop_var81872 = 0; loop_var81872 < 24; loop_var81872++) {
      int mapped_idx;
      mapped_idx = loop_var81872;
      int al_index_name_symbol;
      al_index_name_symbol = d_z5044dinitial_mpg81658[loop_var81872];
      d_z5044dinitial_der81867[mapped_idx] = (((get_dfdx_var(al_index_name_symbol, 24, dx4794dinitial_inv_mpg81653, dx4794dinitial_der81832) * rot81812[2]) + (get_dfdx_var(al_index_name_symbol, 24, dy4843dinitial_inv_mpg81681, dy4843dinitial_der81839) * rot81812[5])) + (get_dfdx_var(al_index_name_symbol, 24, dz4892dinitial_inv_mpg81675, dz4892dinitial_der81846) * rot81812[8]));
    }
    _z81868 = (((rot81812[2] * x81833) + (rot81812[5] * y81840)) + (rot81812[8] * z81847));
    long double _r281875;
    long double d_r25100dinitial_der81874[24] = { 0.0 };
    for (int loop_var81879 = 0; loop_var81879 < 24; loop_var81879++) {
      int mapped_idx;
      mapped_idx = loop_var81879;
      int al_index_name_symbol;
      al_index_name_symbol = d_r25100dinitial_mpg81670[loop_var81879];
      d_r25100dinitial_der81874[mapped_idx] = ((((_x81854 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853)) + (_x81854 * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853))) + ((_y81861 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860)) + (_y81861 * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860)))) + ((_z81868 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867))));
    }
    _r281875 = (((_x81854 * _x81854) + (_y81861 * _y81861)) + (_z81868 * _z81868));
    long double _r81882;
    long double d_r5143dinitial_der81881[24] = { 0.0 };
    for (int loop_var81886 = 0; loop_var81886 < 24; loop_var81886++) {
      int mapped_idx;
      mapped_idx = loop_var81886;
      int al_index_name_symbol;
      al_index_name_symbol = d_r5143dinitial_mpg81690[loop_var81886];
      d_r5143dinitial_der81881[mapped_idx] = (0.5 * (sqrtl((1.0 / _r281875)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874)));
    }
    _r81882 = sqrtl(_r281875);
    long double _r381889;
    long double d_r35167dinitial_der81888[24] = { 0.0 };
    for (int loop_var81893 = 0; loop_var81893 < 24; loop_var81893++) {
      int mapped_idx;
      mapped_idx = loop_var81893;
      int al_index_name_symbol;
      al_index_name_symbol = d_r35167dinitial_mpg81672[loop_var81893];
      d_r35167dinitial_der81888[mapped_idx] = ((_r81882 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874)) + (_r281875 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81691, d_r5143dinitial_der81881)));
    }
    _r381889 = (_r81882 * _r281875);
    long double _r481896;
    long double d_r45191dinitial_der81895[24] = { 0.0 };
    for (int loop_var81900 = 0; loop_var81900 < 24; loop_var81900++) {
      int mapped_idx;
      mapped_idx = loop_var81900;
      int al_index_name_symbol;
      al_index_name_symbol = d_r45191dinitial_mpg81698[loop_var81900];
      d_r45191dinitial_der81895[mapped_idx] = ((_r281875 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874)) + (_r281875 * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874)));
    }
    _r481896 = (_r281875 * _r281875);
    long double _r581903;
    long double d_r55216dinitial_der81902[24] = { 0.0 };
    for (int loop_var81907 = 0; loop_var81907 < 24; loop_var81907++) {
      int mapped_idx;
      mapped_idx = loop_var81907;
      int al_index_name_symbol;
      al_index_name_symbol = d_r55216dinitial_mpg81646[loop_var81907];
      d_r55216dinitial_der81902[mapped_idx] = ((_r481896 * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81691, d_r5143dinitial_der81881)) + (_r81882 * get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895)));
    }
    _r581903 = (_r481896 * _r81882);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81696[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81813[((24 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r81882 * _r81882)) * get_dfdx_var(al_index_name_symbol, 24, d_r5143dinitial_inv_mpg81691, d_r5143dinitial_der81881))));
        }
      }
    }
    v_081814[0] = (1.0e0 / _r81882);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg81696[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der81813[((24 * 1) + mapped_idx)] = ((((get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * 1.7320508075688772936e0) * _r381889) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81673, d_r35167dinitial_der81888) * (_z81868 * 1.7320508075688772936e0))) / (_r381889 * _r381889));
        }
      }
    }
    v_081814[1] = ((1.7320508075688772936e0 * _z81868) / _r381889);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81660[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81815[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853) * _r381889) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81673, d_r35167dinitial_der81888) * -_x81854)) / (_r381889 * _r381889));
        }
      }
    }
    dv_0_dx81816[0] = (-_x81854 / _r381889);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81692[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81817[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860) * _r381889) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81673, d_r35167dinitial_der81888) * -_y81861)) / (_r381889 * _r381889));
        }
      }
    }
    dv_0_dy81818[0] = (-_y81861 / _r381889);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81666[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81819[((24 * 0) + mapped_idx)] = (((-get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * _r381889) - (get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81673, d_r35167dinitial_der81888) * -_z81868)) / (_r381889 * _r381889));
        }
      }
    }
    dv_0_dz81820[0] = (-_z81868 / _r381889);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg81660[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der81815[((24 * 1) + mapped_idx)] = ((((((_z81868 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853)) + (_x81854 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * -5.1961524227066318805e0))) * _r581903) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81647, d_r55216dinitial_der81902) * ((_z81868 * -5.1961524227066318805e0) * _x81854))) / (_r581903 * _r581903));
        }
      }
    }
    dv_0_dx81816[1] = (((-5.1961524227066318805e0 * _z81868) * _x81854) / _r581903);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg81692[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der81817[((24 * 1) + mapped_idx)] = ((((((_z81868 * -5.1961524227066318805e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860)) + (_y81861 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * -5.1961524227066318805e0))) * _r581903) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81647, d_r55216dinitial_der81902) * ((_z81868 * -5.1961524227066318805e0) * _y81861))) / (_r581903 * _r581903));
        }
      }
    }
    dv_0_dy81818[1] = (((-5.1961524227066318805e0 * _z81868) * _y81861) / _r581903);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg81666[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der81819[((24 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r381889 * _r381889)) * get_dfdx_var(al_index_name_symbol, 24, d_r35167dinitial_inv_mpg81673, d_r35167dinitial_der81888)))) - ((((((_z81868 * 3.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * 3.0e0))) * _r581903) - (get_dfdx_var(al_index_name_symbol, 24, d_r55216dinitial_inv_mpg81647, d_r55216dinitial_der81902) * ((_z81868 * 3.0e0) * _z81868))) / (_r581903 * _r581903))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz81820[1] = (1.7320508075688772936e0 * ((1.0e0 / _r381889) - (((3.0e0 * _z81868) * _z81868) / _r581903)));
    for (int n81941 = 2; n81941 < 7; n81941++) {
      long double coef181943;
      coef181943 = sqrtl(((((2.0e0 * ((long double) n81941)) - 1.0e0) * ((long double) ((2 * n81941) + 1))) / ((long double) (n81941 * n81941))));
      long double coef281947;
      coef281947 = sqrtl(((((((long double) n81941) - 1.0e0) * ((long double) (n81941 - 1))) * ((long double) ((2 * n81941) + 1))) / ((long double) ((n81941 * n81941) * ((2 * n81941) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81941));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_04620dinitial_mpg81696[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dv_04620dinitial_der81813[((24 * n81941) + mapped_idx)] = (((((((v_081814[(n81941 - 1)] * coef181943) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * (get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * coef181943))) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * ((v_081814[(n81941 - 1)] * coef181943) * _z81868))) / (_r281875 * _r281875)) - ((((get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * coef281947) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * (v_081814[(n81941 - 2)] * coef281947))) / (_r281875 * _r281875)));
          }
        }
      }
      v_081814[n81941] = ((((coef181943 * v_081814[(n81941 - 1)]) * _z81868) / _r281875) - ((coef281947 * v_081814[(n81941 - 2)]) / _r281875));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81941));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx4638dinitial_mpg81660[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dx4638dinitial_der81815[((24 * n81941) + mapped_idx)] = ((((_z81868 * coef181943) * ((((get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81661, ddv_0_dx4638dinitial_der81815) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * dv_0_dx81816[(n81941 - 1)])) / (_r281875 * _r281875)) - ((((((v_081814[(n81941 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853)) + (_x81854 * (get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((v_081814[(n81941 - 1)] * 2.0e0) * _x81854))) / (_r481896 * _r481896)))) + (((dv_0_dx81816[(n81941 - 1)] / _r281875) - (((v_081814[(n81941 - 1)] * 2.0e0) * _x81854) / _r481896)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * coef181943))) - (((((get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81661, ddv_0_dx4638dinitial_der81815) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * dv_0_dx81816[(n81941 - 2)])) / (_r281875 * _r281875)) - ((((((v_081814[(n81941 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_x4942dinitial_inv_mpg81687, d_x4942dinitial_der81853)) + (_x81854 * (get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((v_081814[(n81941 - 2)] * 2.0e0) * _x81854))) / (_r481896 * _r481896))) * coef281947));
          }
        }
      }
      dv_0_dx81816[n81941] = (((coef181943 * _z81868) * ((dv_0_dx81816[(n81941 - 1)] / _r281875) - (((v_081814[(n81941 - 1)] * 2.0e0) * _x81854) / _r481896))) - (coef281947 * ((dv_0_dx81816[(n81941 - 2)] / _r281875) - (((v_081814[(n81941 - 2)] * 2.0e0) * _x81854) / _r481896))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81941));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy4660dinitial_mpg81692[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dy4660dinitial_der81817[((24 * n81941) + mapped_idx)] = ((((_z81868 * coef181943) * ((((get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81693, ddv_0_dy4660dinitial_der81817) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * dv_0_dy81818[(n81941 - 1)])) / (_r281875 * _r281875)) - ((((((v_081814[(n81941 - 1)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860)) + (_y81861 * (get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((v_081814[(n81941 - 1)] * 2.0e0) * _y81861))) / (_r481896 * _r481896)))) + (((dv_0_dy81818[(n81941 - 1)] / _r281875) - (((v_081814[(n81941 - 1)] * 2.0e0) * _y81861) / _r481896)) * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * coef181943))) - (((((get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81693, ddv_0_dy4660dinitial_der81817) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * dv_0_dy81818[(n81941 - 2)])) / (_r281875 * _r281875)) - ((((((v_081814[(n81941 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_y4993dinitial_inv_mpg81655, d_y4993dinitial_der81860)) + (_y81861 * (get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((v_081814[(n81941 - 2)] * 2.0e0) * _y81861))) / (_r481896 * _r481896))) * coef281947));
          }
        }
      }
      dv_0_dy81818[n81941] = (((coef181943 * _z81868) * ((dv_0_dy81818[(n81941 - 1)] / _r281875) - (((v_081814[(n81941 - 1)] * 2.0e0) * _y81861) / _r481896))) - (coef281947 * ((dv_0_dy81818[(n81941 - 2)] / _r281875) - (((v_081814[(n81941 - 2)] * 2.0e0) * _y81861) / _r481896))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n81941));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz4682dinitial_mpg81666[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dz4682dinitial_der81819[((24 * n81941) + mapped_idx)] = ((((((((dv_0_dz81820[(n81941 - 1)] * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81667, ddv_0_dz4682dinitial_der81819))) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * (dv_0_dz81820[(n81941 - 1)] * _z81868))) / (_r281875 * _r281875)) + ((v_081814[(n81941 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r281875 * _r281875)) * get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874)))) - ((((((_z81868 * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * (get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((_z81868 * 2.0e0) * _z81868))) / (_r481896 * _r481896)))) + (((1.0e0 / _r281875) - (((_z81868 * 2.0e0) * _z81868) / _r481896)) * get_dfdx_cell((n81941 - 1), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813)))) * coef181943) - (((((get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81667, ddv_0_dz4682dinitial_der81819) * _r281875) - (get_dfdx_var(al_index_name_symbol, 24, d_r25100dinitial_inv_mpg81671, d_r25100dinitial_der81874) * dv_0_dz81820[(n81941 - 2)])) / (_r281875 * _r281875)) - ((((((v_081814[(n81941 - 2)] * 2.0e0) * get_dfdx_var(al_index_name_symbol, 24, d_z5044dinitial_inv_mpg81659, d_z5044dinitial_der81867)) + (_z81868 * (get_dfdx_cell((n81941 - 2), 24, al_index_name_symbol, 24, dv_04620dinitial_inv_mpg81697, dv_04620dinitial_der81813) * 2.0e0))) * _r481896) - (get_dfdx_var(al_index_name_symbol, 24, d_r45191dinitial_inv_mpg81699, d_r45191dinitial_der81895) * ((v_081814[(n81941 - 2)] * 2.0e0) * _z81868))) / (_r481896 * _r481896))) * coef281947));
          }
        }
      }
      dv_0_dz81820[n81941] = ((coef181943 * (((dv_0_dz81820[(n81941 - 1)] * _z81868) / _r281875) + (v_081814[(n81941 - 1)] * ((1.0e0 / _r281875) - (((2.0e0 * _z81868) * _z81868) / _r481896))))) - (coef281947 * ((dv_0_dz81820[(n81941 - 2)] / _r281875) - (((v_081814[(n81941 - 2)] * 2.0e0) * _z81868) / _r481896))));
    }
    long double dpotential_dx81968;
    long double ddpotential_dx6387dinitial_der81967[24] = { 0.0 };
    for (int loop_var81972 = 0; loop_var81972 < 24; loop_var81972++) {
      int mapped_idx;
      mapped_idx = loop_var81972;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx6387dinitial_mpg81662[loop_var81972];
      ddpotential_dx6387dinitial_der81967[mapped_idx] = 0.0;
    }
    dpotential_dx81968 = ((long double) 0);
    long double dpotential_dy81974;
    long double ddpotential_dy6414dinitial_der81973[24] = { 0.0 };
    for (int loop_var81978 = 0; loop_var81978 < 24; loop_var81978++) {
      int mapped_idx;
      mapped_idx = loop_var81978;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy6414dinitial_mpg81702[loop_var81978];
      ddpotential_dy6414dinitial_der81973[mapped_idx] = 0.0;
    }
    dpotential_dy81974 = ((long double) 0);
    long double dpotential_dz81980;
    long double ddpotential_dz6441dinitial_der81979[24] = { 0.0 };
    for (int loop_var81984 = 0; loop_var81984 < 24; loop_var81984++) {
      int mapped_idx;
      mapped_idx = loop_var81984;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz6441dinitial_mpg81694[loop_var81984];
      ddpotential_dz6441dinitial_der81979[mapped_idx] = 0.0;
    }
    dpotential_dz81980 = ((long double) 0);
    for (int n81985 = 2; n81985 < 7; n81985++) {
      for (int loop_var81990 = 0; loop_var81990 < 24; loop_var81990++) {
        int mapped_idx;
        mapped_idx = loop_var81990;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx6387dinitial_mpg81662[loop_var81990];
        ddpotential_dx6387dinitial_der81967[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81663, ddpotential_dx6387dinitial_der81967) + (get_dfdx_cell(n81985, 24, al_index_name_symbol, 24, ddv_0_dx4638dinitial_inv_mpg81661, ddv_0_dx4638dinitial_der81815) * central_grav81573[n81985]));
      }
      dpotential_dx81968 = (dpotential_dx81968 + (dv_0_dx81816[n81985] * central_grav81573[n81985]));
      for (int loop_var81995 = 0; loop_var81995 < 24; loop_var81995++) {
        int mapped_idx;
        mapped_idx = loop_var81995;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy6414dinitial_mpg81702[loop_var81995];
        ddpotential_dy6414dinitial_der81973[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81703, ddpotential_dy6414dinitial_der81973) + (get_dfdx_cell(n81985, 24, al_index_name_symbol, 24, ddv_0_dy4660dinitial_inv_mpg81693, ddv_0_dy4660dinitial_der81817) * central_grav81573[n81985]));
      }
      dpotential_dy81974 = (dpotential_dy81974 + (dv_0_dy81818[n81985] * central_grav81573[n81985]));
      for (int loop_var82000 = 0; loop_var82000 < 24; loop_var82000++) {
        int mapped_idx;
        mapped_idx = loop_var82000;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz6441dinitial_mpg81694[loop_var82000];
        ddpotential_dz6441dinitial_der81979[mapped_idx] = (get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81695, ddpotential_dz6441dinitial_der81979) + (get_dfdx_cell(n81985, 24, al_index_name_symbol, 24, ddv_0_dz4682dinitial_inv_mpg81667, ddv_0_dz4682dinitial_der81819) * central_grav81573[n81985]));
      }
      dpotential_dz81980 = (dpotential_dz81980 + (dv_0_dz81820[n81985] * central_grav81573[n81985]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81682[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81806[((24 * 0) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dx6387dinitial_inv_mpg81663, ddpotential_dx6387dinitial_der81967);
        }
      }
    }
    grav_acc_local81807[0] = dpotential_dx81968;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81682[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81806[((24 * 1) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dy6414dinitial_inv_mpg81703, ddpotential_dy6414dinitial_der81973);
        }
      }
    }
    grav_acc_local81807[1] = dpotential_dy81974;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 2));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg81682[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der81806[((24 * 2) + mapped_idx)] = get_dfdx_var(al_index_name_symbol, 24, ddpotential_dz6441dinitial_inv_mpg81695, ddpotential_dz6441dinitial_der81979);
        }
      }
    }
    grav_acc_local81807[2] = dpotential_dz81980;
    for (int k82014 = 0; k82014 < 3; k82014++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * k82014));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc4569dinitial_mpg81688[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dgrav_acc4569dinitial_der81808[((24 * k82014) + mapped_idx)] = ((((get_dfdx_cell(0, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81683, dgrav_acc_local4544dinitial_der81806) * rot81812[((3 * k82014) + 0)]) + (get_dfdx_cell(1, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81683, dgrav_acc_local4544dinitial_der81806) * rot81812[((3 * k82014) + 1)])) + (get_dfdx_cell(2, 24, al_index_name_symbol, 24, dgrav_acc_local4544dinitial_inv_mpg81683, dgrav_acc_local4544dinitial_der81806) * rot81812[((3 * k82014) + 2)])) / 2.2838315556293922983e-07);
          }
        }
      }
      grav_acc81809[k82014] = ((((rot81812[((3 * k82014) + 0)] * grav_acc_local81807[0]) + (rot81812[((3 * k82014) + 1)] * grav_acc_local81807[1])) + (rot81812[((3 * k82014) + 2)] * grav_acc_local81807[2])) / 2.2838315556293922983e-07);
    }
    for (int slice_idx = 0; slice_idx < (((i81830 * 3) + 3) - (i81830 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i81830 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg81700[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der81704[((24 * ((i81830 * 3) + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i81830 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81701, dsat_acc2828dinitial_der81704) + (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81689, dgrav_acc4569dinitial_der81808) * 2.8247609439046209905e-07));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i81830 * 3) + 3) - (i81830 * 3)); slice_idx++) {
      sat_acc81705[((i81830 * 3) + slice_idx)] = (sat_acc81705[((i81830 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc81809[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg81668[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der81706[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81669, dcentral_acc2846dinitial_der81706) - (get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dgrav_acc4569dinitial_inv_mpg81689, dgrav_acc4569dinitial_der81808) * sat_gms81644[i81830]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc81707[(0 + slice_idx)] = (central_acc81707[(0 + slice_idx)] - (sat_gms81644[i81830] * grav_acc81809[(0 + slice_idx)]));
    }
  }
  for (int i82028 = 0; i82028 < 4; i82028++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc4588dinitial_mpg81648[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dacc4588dinitial_der81810[((24 * (0 + slice_idx)) + mapped_idx)] = (get_dfdx_cell(((i82028 * 3) + slice_idx), 24, al_index_name_symbol, 24, dsat_acc2828dinitial_inv_mpg81701, dsat_acc2828dinitial_der81704) - get_dfdx_cell((0 + slice_idx), 24, al_index_name_symbol, 24, dcentral_acc2846dinitial_inv_mpg81669, dcentral_acc2846dinitial_der81706));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc81811[(0 + slice_idx)] = (sat_acc81705[((i82028 * 3) + slice_idx)] - central_acc81707[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82028 * 6) + 3) - (i82028 * 6)); slice_idx++) {
        jupsatsystem81638[((i82028 * 6) + slice_idx)] = state81710[(((i82028 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82028 * 6) + 6) - ((i82028 * 6) + 3)); slice_idx++) {
        jupsatsystem81638[(((i82028 * 6) + 3) + slice_idx)] = acc81811[(0 + slice_idx)];
      }
    }
    for (int j82040 = 0; j82040 < 3; j82040++) {
      {
        for (int slice_idx = 0; slice_idx < ((24 + (((((i82028 * 6) + j82040) + 1) * 4) * 6)) - (24 + ((((i82028 * 6) + j82040) * 4) * 6))); slice_idx++) {
          jupsatsystem81638[((24 + ((((i82028 * 6) + j82040) * 4) * 6)) + slice_idx)] = state_and_derivatives81640[((24 + (((((i82028 * 6) + j82040) + 3) * 4) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((24 + (((((i82028 * 6) + j82040) + 4) * 4) * 6)) - (24 + (((((i82028 * 6) + j82040) + 3) * 4) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (24 + (((((i82028 * 6) + j82040) + 3) * 4) * 6)));
        jupsatsystem81638[func_slice_idx] = get_dfdx_cell(j82040, 24, al_index_name_symbol, 24, dacc4588dinitial_inv_mpg81649, dacc4588dinitial_der81810);
      }
    }
  }
  return 0;
}

int jupsatsystem_noderiv(long double *restrict jupsatsystem_noderiv82047, long double t82048, long double *restrict state82049, long double *restrict central_pos82050, long double *restrict perturb_gms82051, long double *restrict perturb_pos82052, long double *restrict sat_gms82053) {
  long double sat_acc82055[12] = { 0.0 };
  long double central_acc82056[3] = { 0.0 };
  long double dist282057;
  long double dist382058;
  for (int i82059 = 0; i82059 < 4; i82059++) {
    long double r82061[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82061[(0 + slice_idx)] = state82049[((i82059 * 6) + slice_idx)];
      }
    }
    dist282057 = (((r82061[0] * r82061[0]) + (r82061[1] * r82061[1])) + (r82061[2] * r82061[2]));
    dist382058 = (dist282057 * sqrtl(dist282057));
    {
      for (int slice_idx = 0; slice_idx < (((i82059 * 3) + 3) - (i82059 * 3)); slice_idx++) {
        sat_acc82055[((i82059 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r82061[(0 + slice_idx)]) / dist382058);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82056[(0 + slice_idx)] = (central_acc82056[(0 + slice_idx)] + ((sat_gms82053[i82059] * r82061[(0 + slice_idx)]) / dist382058));
      }
    }
  }
  for (int i82077 = 0; i82077 < 4; i82077++) {
    long double r82079[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r82079[(0 + slice_idx)] = (perturb_pos82052[((i82077 * 3) + slice_idx)] - central_pos82050[(0 + slice_idx)]);
      }
    }
    dist282057 = (((r82079[0] * r82079[0]) + (r82079[1] * r82079[1])) + (r82079[2] * r82079[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82056[(0 + slice_idx)] = (central_acc82056[(0 + slice_idx)] + (((perturb_gms82051[i82077] * r82079[(0 + slice_idx)]) / dist282057) / sqrtl(dist282057)));
      }
    }
    for (int j82089 = 0; j82089 < 4; j82089++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82079[(0 + slice_idx)] = (perturb_pos82052[((i82077 * 3) + slice_idx)] - (state82049[((j82089 * 6) + slice_idx)] + central_pos82050[(0 + slice_idx)]));
        }
      }
      dist282057 = (((r82079[0] * r82079[0]) + (r82079[1] * r82079[1])) + (r82079[2] * r82079[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j82089 * 3) + 3) - (j82089 * 3)); slice_idx++) {
          sat_acc82055[((j82089 * 3) + slice_idx)] = (sat_acc82055[((j82089 * 3) + slice_idx)] + (((perturb_gms82051[i82077] * r82079[(0 + slice_idx)]) / dist282057) / sqrtl(dist282057)));
        }
      }
    }
  }
  for (int i82100 = 1; i82100 < 4; i82100++) {
    long double r82102[3] = { 0.0 };
    for (int j82103 = 0; j82103 < i82100; j82103++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r82102[(0 + slice_idx)] = (state82049[((j82103 * 6) + slice_idx)] - state82049[((i82100 * 6) + slice_idx)]);
        }
      }
      dist282057 = (((r82102[0] * r82102[0]) + (r82102[1] * r82102[1])) + (r82102[2] * r82102[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i82100 * 3) + 3) - (i82100 * 3)); slice_idx++) {
          sat_acc82055[((i82100 * 3) + slice_idx)] = (sat_acc82055[((i82100 * 3) + slice_idx)] + (((sat_gms82053[j82103] * r82102[(0 + slice_idx)]) / dist282057) / sqrtl(dist282057)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j82103 * 3) + 3) - (j82103 * 3)); slice_idx++) {
          sat_acc82055[((j82103 * 3) + slice_idx)] = (sat_acc82055[((j82103 * 3) + slice_idx)] - (((sat_gms82053[i82100] * r82102[(0 + slice_idx)]) / dist282057) / sqrtl(dist282057)));
        }
      }
    }
  }
  long double grav_acc_local82117[3] = { 0.0 };
  long double grav_acc82118[3] = { 0.0 };
  long double acc82119[3] = { 0.0 };
  long double rot82120[9] = { 0.0 };
  long double v_082121[7] = { 0.0 };
  long double dv_0_dx82122[7] = { 0.0 };
  long double dv_0_dy82123[7] = { 0.0 };
  long double dv_0_dz82124[7] = { 0.0 };
  {
    long double jupiter_rotation_matrix82133[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix82133, t82048);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot82120[(0 + slice_idx)] = jupiter_rotation_matrix82133[slice_idx];
    }
  }
  for (int i82134 = 0; i82134 < 4; i82134++) {
    long double x82136;
    x82136 = (state82049[((i82134 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y82140;
    y82140 = (state82049[((i82134 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z82144;
    z82144 = (state82049[((i82134 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x82148;
    _x82148 = (((rot82120[0] * x82136) + (rot82120[3] * y82140)) + (rot82120[6] * z82144));
    long double _y82152;
    _y82152 = (((rot82120[1] * x82136) + (rot82120[4] * y82140)) + (rot82120[7] * z82144));
    long double _z82156;
    _z82156 = (((rot82120[2] * x82136) + (rot82120[5] * y82140)) + (rot82120[8] * z82144));
    long double _r282160;
    _r282160 = (((_x82148 * _x82148) + (_y82152 * _y82152)) + (_z82156 * _z82156));
    long double r82164;
    r82164 = sqrtl(_r282160);
    long double _r382168;
    _r382168 = (r82164 * _r282160);
    long double _r482172;
    _r482172 = (_r282160 * _r282160);
    long double _r582176;
    _r582176 = (_r482172 * r82164);
    v_082121[0] = (1.0e0 / r82164);
    v_082121[1] = ((1.7320508075688772936e0 * _z82156) / _r382168);
    dv_0_dx82122[0] = (-_x82148 / _r382168);
    dv_0_dy82123[0] = (-_y82152 / _r382168);
    dv_0_dz82124[0] = (-_z82156 / _r382168);
    dv_0_dx82122[1] = (((-5.1961524227066318805e0 * _z82156) * _x82148) / _r582176);
    dv_0_dy82123[1] = (((-5.1961524227066318805e0 * _z82156) * _y82152) / _r582176);
    dv_0_dz82124[1] = (1.7320508075688772936e0 * ((1.0e0 / _r382168) - (((3.0e0 * _z82156) * _z82156) / _r582176)));
    for (int n82204 = 2; n82204 < 7; n82204++) {
      long double coef182206;
      coef182206 = sqrtl(((((2.0e0 * ((long double) n82204)) - 1.0e0) * ((long double) ((2 * n82204) + 1))) / ((long double) (n82204 * n82204))));
      long double coef282210;
      coef282210 = sqrtl(((((((long double) n82204) - 1.0e0) * ((long double) (n82204 - 1))) * ((long double) ((2 * n82204) + 1))) / ((long double) ((n82204 * n82204) * ((2 * n82204) - 3)))));
      v_082121[n82204] = ((((coef182206 * v_082121[(n82204 - 1)]) * _z82156) / _r282160) - ((coef282210 * v_082121[(n82204 - 2)]) / _r282160));
      dv_0_dx82122[n82204] = (((coef182206 * _z82156) * ((dv_0_dx82122[(n82204 - 1)] / _r282160) - (((v_082121[(n82204 - 1)] * 2.0e0) * _x82148) / _r482172))) - (coef282210 * ((dv_0_dx82122[(n82204 - 2)] / _r282160) - (((v_082121[(n82204 - 2)] * 2.0e0) * _x82148) / _r482172))));
      dv_0_dy82123[n82204] = (((coef182206 * _z82156) * ((dv_0_dy82123[(n82204 - 1)] / _r282160) - (((v_082121[(n82204 - 1)] * 2.0e0) * _y82152) / _r482172))) - (coef282210 * ((dv_0_dy82123[(n82204 - 2)] / _r282160) - (((v_082121[(n82204 - 2)] * 2.0e0) * _y82152) / _r482172))));
      dv_0_dz82124[n82204] = ((coef182206 * (((dv_0_dz82124[(n82204 - 1)] * _z82156) / _r282160) + (v_082121[(n82204 - 1)] * ((1.0e0 / _r282160) - (((2.0e0 * _z82156) * _z82156) / _r482172))))) - (coef282210 * ((dv_0_dz82124[(n82204 - 2)] / _r282160) - (((v_082121[(n82204 - 2)] * 2.0e0) * _z82156) / _r482172))));
    }
    long double dpotential_dx82226;
    dpotential_dx82226 = ((long double) 0);
    long double dpotential_dy82230;
    dpotential_dy82230 = ((long double) 0);
    long double dpotential_dz82234;
    dpotential_dz82234 = ((long double) 0);
    for (int n82238 = 2; n82238 < 7; n82238++) {
      dpotential_dx82226 = (dpotential_dx82226 + (dv_0_dx82122[n82238] * central_grav81573[n82238]));
      dpotential_dy82230 = (dpotential_dy82230 + (dv_0_dy82123[n82238] * central_grav81573[n82238]));
      dpotential_dz82234 = (dpotential_dz82234 + (dv_0_dz82124[n82238] * central_grav81573[n82238]));
    }
    grav_acc_local82117[0] = dpotential_dx82226;
    grav_acc_local82117[1] = dpotential_dy82230;
    grav_acc_local82117[2] = dpotential_dz82234;
    for (int k82258 = 0; k82258 < 3; k82258++) {
      grav_acc82118[k82258] = ((((rot82120[((3 * k82258) + 0)] * grav_acc_local82117[0]) + (rot82120[((3 * k82258) + 1)] * grav_acc_local82117[1])) + (rot82120[((3 * k82258) + 2)] * grav_acc_local82117[2])) / 2.2838315556293922983e-07);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82134 * 3) + 3) - (i82134 * 3)); slice_idx++) {
        sat_acc82055[((i82134 * 3) + slice_idx)] = (sat_acc82055[((i82134 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc82118[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc82056[(0 + slice_idx)] = (central_acc82056[(0 + slice_idx)] - (sat_gms82053[i82134] * grav_acc82118[(0 + slice_idx)]));
      }
    }
  }
  for (int i82269 = 0; i82269 < 4; i82269++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc82119[(0 + slice_idx)] = (sat_acc82055[((i82269 * 3) + slice_idx)] - central_acc82056[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82269 * 6) + 3) - (i82269 * 6)); slice_idx++) {
        jupsatsystem_noderiv82047[((i82269 * 6) + slice_idx)] = state82049[(((i82269 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i82269 * 6) + 6) - ((i82269 * 6) + 3)); slice_idx++) {
        jupsatsystem_noderiv82047[(((i82269 * 6) + 3) + slice_idx)] = acc82119[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

