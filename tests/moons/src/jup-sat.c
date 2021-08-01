#include <math.h>

const long double pole_ra77968[2] = { 268.05659500000001572e0, -0.006498999999999999569e0 };
const long double pole_dec77969[2] = { 64.49530300000000693e0, 0.0024130000000000002142e0 };
const long double pm77970[2] = { 284.94999999999998863e0, 870.5359999999999445e0 };
const long double nut_prec_ra77971[5] = { 0.00011699999999999999788e0, 0.00093800000000000003167e0, 0.0014319999999999998962e0, 3.000000000000000076e-05, 0.0021500000000000000014e0 };
const long double nut_prec_dec77972[5] = { 5.0000000000000002396e-05, 0.00040400000000000000798e0, 0.0006170000000000000354e0, -1.29999999999999992e-05, 0.0009259999999999999568e0 };
const long double Jabcde_077973[5] = { 99.36071400000000153e0, 175.89536899999998809e0, 300.32316200000002482e0, 114.01230499999999779e0, 49.511251000000001454e0 };
const long double Jabcde_T77974[5] = { 4850.4045999999998457e0, 1191.9604999999999109e0, 262.54750000000001364e0, 6070.247599999999693e0, 64.29999999999999716e0 };
const long double central_grav77975[7] = { 0.0e0, 0.0e0, -0.0065724808672554692115e0, 0.0e0, 0.00019554099999999999462e0, 0.0e0, -9.4975767597683728686e-06 };

int jupiter_rotation_matrix(long double *restrict jupiter_rotation_matrix77976, long double t77977) {
  long double T77979;
  T77979 = (t77977 / 36525.0e0);
  long double alpha_077983;
  alpha_077983 = (268.05659500000001572e0 + (-0.006498999999999999569e0 * T77979));
  long double delta_077987;
  delta_077987 = (64.49530300000000693e0 + (0.0024130000000000002142e0 * T77979));
  long double W77991;
  W77991 = (((284.94999999999998863e0 + (870.5359999999999445e0 * t77977)) * 3.141592653589793116e0) / 180.0e0);
  for (int i77995 = 0; i77995 < 5; i77995++) {
    long double J77997;
    J77997 = (((Jabcde_077973[i77995] + (Jabcde_T77974[i77995] * T77979)) * 3.141592653589793116e0) / 180.0e0);
    alpha_077983 = (alpha_077983 + (nut_prec_ra77971[i77995] * sinl(J77997)));
    delta_077987 = (delta_077987 + (nut_prec_dec77972[i77995] * cosl(J77997)));
  }
  alpha_077983 = (alpha_077983 * 0.017453292519943295088e0);
  delta_077987 = (delta_077987 * 0.017453292519943295088e0);
  jupiter_rotation_matrix77976[0] = ((-sinl(alpha_077983) * cosl(W77991)) - ((cosl(alpha_077983) * sinl(delta_077987)) * sinl(W77991)));
  jupiter_rotation_matrix77976[1] = ((sinl(alpha_077983) * sinl(W77991)) - ((cosl(alpha_077983) * sinl(delta_077987)) * cosl(W77991)));
  jupiter_rotation_matrix77976[2] = (cosl(alpha_077983) * cosl(delta_077987));
  jupiter_rotation_matrix77976[3] = ((cosl(alpha_077983) * cosl(W77991)) - ((sinl(alpha_077983) * sinl(delta_077987)) * sinl(W77991)));
  jupiter_rotation_matrix77976[4] = ((-cosl(alpha_077983) * sinl(W77991)) - ((sinl(alpha_077983) * sinl(delta_077987)) * cosl(W77991)));
  jupiter_rotation_matrix77976[5] = (cosl(delta_077987) * sinl(alpha_077983));
  jupiter_rotation_matrix77976[6] = (cosl(delta_077987) * sinl(W77991));
  jupiter_rotation_matrix77976[7] = (cosl(delta_077987) * cosl(W77991));
  jupiter_rotation_matrix77976[8] = sinl(delta_077987);
  return 0;
}

int jupsatsystem(long double *restrict jupsatsystem78040, long double t78041, long double *restrict state_and_derivatives78042, long double *restrict central_pos78043, long double *restrict perturb_gms78044, long double *restrict perturb_pos78045, long double *restrict sat_gms78046) {
  static int ddist33107dinitial_mpg78048[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist33107dinitial_inv_mpg78049[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_mpg78050[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dx6387dinitial_inv_mpg78051[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dstate2924dinitial_mpg78052[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dstate2924dinitial_inv_mpg78053[576] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dcentral_acc2846dinitial_mpg78054[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dcentral_acc2846dinitial_inv_mpg78055[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_mpg78056[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_z5044dinitial_inv_mpg78057[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_mpg78058[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dx4638dinitial_inv_mpg78059[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_mpg78060[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r5143dinitial_inv_mpg78061[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_mpg78062[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dsat_acc2828dinitial_inv_mpg78063[288] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_mpg78064[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dz4682dinitial_inv_mpg78065[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_mpg78066[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr4219dinitial_inv_mpg78067[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_mpg78068[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dz6441dinitial_inv_mpg78069[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_mpg78070[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r35167dinitial_inv_mpg78071[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dv_04620dinitial_mpg78072[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dv_04620dinitial_inv_mpg78073[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  static int dgrav_acc_local4544dinitial_mpg78074[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc_local4544dinitial_inv_mpg78075[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_mpg78076[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_x4942dinitial_inv_mpg78077[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_mpg78078[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dgrav_acc4569dinitial_inv_mpg78079[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_mpg78080[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddv_0_dy4660dinitial_inv_mpg78081[168] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_mpg78082[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dx4794dinitial_inv_mpg78083[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_mpg78084[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3708dinitial_inv_mpg78085[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_mpg78086[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dr3376dinitial_inv_mpg78087[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_mpg78088[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddpotential_dy6414dinitial_inv_mpg78089[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_mpg78090[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r55216dinitial_inv_mpg78091[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_mpg78092[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r45191dinitial_inv_mpg78093[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_mpg78094[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int ddist23094dinitial_inv_mpg78095[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_mpg78096[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dz4892dinitial_inv_mpg78097[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_mpg78098[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dy4843dinitial_inv_mpg78099[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_mpg78100[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int dacc4588dinitial_inv_mpg78101[72] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_mpg78102[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_y4993dinitial_inv_mpg78103[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_mpg78104[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  static int d_r25100dinitial_inv_mpg78105[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
  long double sat_acc78107[12] = { 0.0 };
  long double dsat_acc2828dinitial_der78106[288] = { 0.0 };
  long double central_acc78109[3] = { 0.0 };
  long double dcentral_acc2846dinitial_der78108[72] = { 0.0 };
  long double state_derivatives_initial78110[576] = { 0.0 };
  long double state78112[24] = { 0.0 };
  long double dstate2924dinitial_der78111[576] = { 0.0 };
  {
    for (int slice_idx = 0; slice_idx < 576; slice_idx++) {
      state_derivatives_initial78110[(0 + slice_idx)] = state_and_derivatives78042[(24 + slice_idx)];
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
        al_index_name_symbol = dstate2924dinitial_mpg78052[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dstate2924dinitial_der78111[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
        }
      }
    }
  }
  for (int slice_idx = 0; slice_idx < 24; slice_idx++) {
    state78112[(0 + slice_idx)] = state_and_derivatives78042[(0 + slice_idx)];
  }
  long double dist278120;
  long double ddist23094dinitial_der78119[24] = { 0.0 };
  long double dist378122;
  long double ddist33107dinitial_der78121[24] = { 0.0 };
  for (int i78123 = 0; i78123 < 24; i78123++) {
    for (int j78125 = 0; j78125 < 24; j78125++) {
      {
        int mapped_idx;
        mapped_idx = ((j78125 < 24) ? dstate2924dinitial_inv_mpg78053[((i78123 * 24) + j78125)] : -1);
        if ((mapped_idx >= 0)) {
          dstate2924dinitial_der78111[((i78123 * 24) + mapped_idx)] = state_derivatives_initial78110[(((i78123 * 4) * 6) + j78125)];
        }
      }
    }
  }
  for (int i78128 = 0; i78128 < 4; i78128++) {
    long double r78131[3] = { 0.0 };
    long double dr3376dinitial_der78130[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3376dinitial_mpg78086[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3376dinitial_der78130[((24 * (0 + slice_idx)) + mapped_idx)] = ((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((i78128 * 6) + slice_idx)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((i78128 * 6) + slice_idx)) + dstate2924dinitial_inv_mpg78053[((24 * ((i78128 * 6) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0);
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r78131[(0 + slice_idx)] = state78112[((i78128 * 6) + slice_idx)];
    }
    for (int loop_var78139 = 0; loop_var78139 < 24; loop_var78139++) {
      int mapped_idx;
      mapped_idx = loop_var78139;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg78094[loop_var78139];
      ddist23094dinitial_der78119[mapped_idx] = ((((r78131[0] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 0) + dr3376dinitial_inv_mpg78087[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78131[0] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 0) + dr3376dinitial_inv_mpg78087[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0))) + ((r78131[1] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 1) + dr3376dinitial_inv_mpg78087[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78131[1] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 1) + dr3376dinitial_inv_mpg78087[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)))) + ((r78131[2] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 2) + dr3376dinitial_inv_mpg78087[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78131[2] * ((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * 2) + dr3376dinitial_inv_mpg78087[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0))));
    }
    dist278120 = (((r78131[0] * r78131[0]) + (r78131[1] * r78131[1])) + (r78131[2] * r78131[2]));
    for (int loop_var78144 = 0; loop_var78144 < 24; loop_var78144++) {
      int mapped_idx;
      mapped_idx = loop_var78144;
      int al_index_name_symbol;
      al_index_name_symbol = ddist33107dinitial_mpg78048[loop_var78144];
      ddist33107dinitial_der78121[mapped_idx] = ((dist278120 * (0.5 * (sqrtl((1.0 / dist278120)) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0)))) + (sqrtl(dist278120) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0)));
    }
    dist378122 = (dist278120 * sqrtl(dist278120));
    for (int slice_idx = 0; slice_idx < (((i78128 * 3) + 3) - (i78128 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i78128 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg78062[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der78106[((24 * ((i78128 * 3) + slice_idx)) + mapped_idx)] = ((((((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * (0 + slice_idx)) + dr3376dinitial_inv_mpg78087[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * -2.8247609439046209905e-07) * dist378122) - (((((al_index_name_symbol < 24) ? ddist33107dinitial_inv_mpg78049[al_index_name_symbol] : -1) >= 0) ? ddist33107dinitial_der78121[((al_index_name_symbol < 24) ? ddist33107dinitial_inv_mpg78049[al_index_name_symbol] : -1)] : 0.0) * (r78131[(0 + slice_idx)] * -2.8247609439046209905e-07))) / (dist378122 * dist378122));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i78128 * 3) + 3) - (i78128 * 3)); slice_idx++) {
      sat_acc78107[((i78128 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r78131[(0 + slice_idx)]) / dist378122);
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg78054[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) + ((((((al_index_name_symbol < 24) ? ((dr3376dinitial_inv_mpg78087[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr3376dinitial_der78130[((24 * (0 + slice_idx)) + dr3376dinitial_inv_mpg78087[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * sat_gms78046[i78128]) * dist378122) - (((((al_index_name_symbol < 24) ? ddist33107dinitial_inv_mpg78049[al_index_name_symbol] : -1) >= 0) ? ddist33107dinitial_der78121[((al_index_name_symbol < 24) ? ddist33107dinitial_inv_mpg78049[al_index_name_symbol] : -1)] : 0.0) * (r78131[(0 + slice_idx)] * sat_gms78046[i78128]))) / (dist378122 * dist378122)));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc78109[(0 + slice_idx)] = (central_acc78109[(0 + slice_idx)] + ((sat_gms78046[i78128] * r78131[(0 + slice_idx)]) / dist378122));
    }
  }
  for (int i78154 = 0; i78154 < 4; i78154++) {
    long double r78157[3] = { 0.0 };
    long double dr3708dinitial_der78156[72] = { 0.0 };
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dr3708dinitial_mpg78084[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dr3708dinitial_der78156[((24 * (0 + slice_idx)) + mapped_idx)] = 0.0;
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      r78157[(0 + slice_idx)] = (perturb_pos78045[((i78154 * 3) + slice_idx)] - central_pos78043[(0 + slice_idx)]);
    }
    for (int loop_var78164 = 0; loop_var78164 < 24; loop_var78164++) {
      int mapped_idx;
      mapped_idx = loop_var78164;
      int al_index_name_symbol;
      al_index_name_symbol = ddist23094dinitial_mpg78094[loop_var78164];
      ddist23094dinitial_der78119[mapped_idx] = ((((r78157[0] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 0) + dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[0] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 0) + dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0))) + ((r78157[1] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 1) + dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[1] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 1) + dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)))) + ((r78157[2] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 2) + dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[2] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 2) + dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0))));
    }
    dist278120 = (((r78157[0] * r78157[0]) + (r78157[1] * r78157[1])) + (r78157[2] * r78157[2]));
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg78054[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) + (((((((((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * (0 + slice_idx)) + dr3708dinitial_inv_mpg78085[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * perturb_gms78044[i78154]) * dist278120) - (((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0) * (r78157[(0 + slice_idx)] * perturb_gms78044[i78154]))) / (dist278120 * dist278120)) * sqrtl(dist278120)) - ((0.5 * (sqrtl((1.0 / dist278120)) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0))) * ((r78157[(0 + slice_idx)] * perturb_gms78044[i78154]) / dist278120))) / (sqrtl(dist278120) * sqrtl(dist278120))));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc78109[(0 + slice_idx)] = (central_acc78109[(0 + slice_idx)] + (((perturb_gms78044[i78154] * r78157[(0 + slice_idx)]) / dist278120) / sqrtl(dist278120)));
    }
    for (int j78170 = 0; j78170 < 4; j78170++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr3708dinitial_mpg78084[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr3708dinitial_der78156[((24 * (0 + slice_idx)) + mapped_idx)] = -((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((j78170 * 6) + slice_idx)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((j78170 * 6) + slice_idx)) + dstate2924dinitial_inv_mpg78053[((24 * ((j78170 * 6) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0);
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r78157[(0 + slice_idx)] = (perturb_pos78045[((i78154 * 3) + slice_idx)] - (state78112[((j78170 * 6) + slice_idx)] + central_pos78043[(0 + slice_idx)]));
      }
      for (int loop_var78179 = 0; loop_var78179 < 24; loop_var78179++) {
        int mapped_idx;
        mapped_idx = loop_var78179;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg78094[loop_var78179];
        ddist23094dinitial_der78119[mapped_idx] = ((((r78157[0] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 0) + dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[0] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 0) + dr3708dinitial_inv_mpg78085[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0))) + ((r78157[1] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 1) + dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[1] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 1) + dr3708dinitial_inv_mpg78085[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)))) + ((r78157[2] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 2) + dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78157[2] * ((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * 2) + dr3708dinitial_inv_mpg78085[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0))));
      }
      dist278120 = (((r78157[0] * r78157[0]) + (r78157[1] * r78157[1])) + (r78157[2] * r78157[2]));
      for (int slice_idx = 0; slice_idx < (((j78170 * 3) + 3) - (j78170 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j78170 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg78062[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der78106[((24 * ((j78170 * 3) + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dsat_acc2828dinitial_inv_mpg78063[((24 * ((j78170 * 3) + slice_idx)) + al_index_name_symbol)] >= 0) ? dsat_acc2828dinitial_der78106[((24 * ((j78170 * 3) + slice_idx)) + dsat_acc2828dinitial_inv_mpg78063[((24 * ((j78170 * 3) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) + (((((((((al_index_name_symbol < 24) ? ((dr3708dinitial_inv_mpg78085[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr3708dinitial_der78156[((24 * (0 + slice_idx)) + dr3708dinitial_inv_mpg78085[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * perturb_gms78044[i78154]) * dist278120) - (((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0) * (r78157[(0 + slice_idx)] * perturb_gms78044[i78154]))) / (dist278120 * dist278120)) * sqrtl(dist278120)) - ((0.5 * (sqrtl((1.0 / dist278120)) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0))) * ((r78157[(0 + slice_idx)] * perturb_gms78044[i78154]) / dist278120))) / (sqrtl(dist278120) * sqrtl(dist278120))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j78170 * 3) + 3) - (j78170 * 3)); slice_idx++) {
        sat_acc78107[((j78170 * 3) + slice_idx)] = (sat_acc78107[((j78170 * 3) + slice_idx)] + (((perturb_gms78044[i78154] * r78157[(0 + slice_idx)]) / dist278120) / sqrtl(dist278120)));
      }
    }
  }
  for (int i78185 = 1; i78185 < 4; i78185++) {
    long double r78188[3] = { 0.0 };
    long double dr4219dinitial_der78187[72] = { 0.0 };
    for (int j78189 = 0; j78189 < i78185; j78189++) {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
          if ((mappings_full_idx_symbol >= 72)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dr4219dinitial_mpg78066[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dr4219dinitial_der78187[((24 * (0 + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((j78189 * 6) + slice_idx)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((j78189 * 6) + slice_idx)) + dstate2924dinitial_inv_mpg78053[((24 * ((j78189 * 6) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) - ((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((i78185 * 6) + slice_idx)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((i78185 * 6) + slice_idx)) + dstate2924dinitial_inv_mpg78053[((24 * ((i78185 * 6) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r78188[(0 + slice_idx)] = (state78112[((j78189 * 6) + slice_idx)] - state78112[((i78185 * 6) + slice_idx)]);
      }
      for (int loop_var78198 = 0; loop_var78198 < 24; loop_var78198++) {
        int mapped_idx;
        mapped_idx = loop_var78198;
        int al_index_name_symbol;
        al_index_name_symbol = ddist23094dinitial_mpg78094[loop_var78198];
        ddist23094dinitial_der78119[mapped_idx] = ((((r78188[0] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 0) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 0) + dr4219dinitial_inv_mpg78067[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78188[0] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 0) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 0) + dr4219dinitial_inv_mpg78067[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0))) + ((r78188[1] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 1) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 1) + dr4219dinitial_inv_mpg78067[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78188[1] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 1) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 1) + dr4219dinitial_inv_mpg78067[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0)))) + ((r78188[2] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 2) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 2) + dr4219dinitial_inv_mpg78067[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0)) + (r78188[2] * ((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * 2) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * 2) + dr4219dinitial_inv_mpg78067[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0))));
      }
      dist278120 = (((r78188[0] * r78188[0]) + (r78188[1] * r78188[1])) + (r78188[2] * r78188[2]));
      for (int slice_idx = 0; slice_idx < (((i78185 * 3) + 3) - (i78185 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((i78185 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg78062[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der78106[((24 * ((i78185 * 3) + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78185 * 3) + slice_idx)) + al_index_name_symbol)] >= 0) ? dsat_acc2828dinitial_der78106[((24 * ((i78185 * 3) + slice_idx)) + dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78185 * 3) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) + (((((((((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * (0 + slice_idx)) + dr4219dinitial_inv_mpg78067[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * sat_gms78046[j78189]) * dist278120) - (((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0) * (r78188[(0 + slice_idx)] * sat_gms78046[j78189]))) / (dist278120 * dist278120)) * sqrtl(dist278120)) - ((0.5 * (sqrtl((1.0 / dist278120)) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0))) * ((r78188[(0 + slice_idx)] * sat_gms78046[j78189]) / dist278120))) / (sqrtl(dist278120) * sqrtl(dist278120))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((i78185 * 3) + 3) - (i78185 * 3)); slice_idx++) {
        sat_acc78107[((i78185 * 3) + slice_idx)] = (sat_acc78107[((i78185 * 3) + slice_idx)] + (((sat_gms78046[j78189] * r78188[(0 + slice_idx)]) / dist278120) / sqrtl(dist278120)));
      }
      for (int slice_idx = 0; slice_idx < (((j78189 * 3) + 3) - (j78189 * 3)); slice_idx++) {
        for (int mapped_idx = 0;; mapped_idx++) {
          int mappings_full_idx_symbol;
          mappings_full_idx_symbol = (mapped_idx + (24 * ((j78189 * 3) + slice_idx)));
          if ((mappings_full_idx_symbol >= 288)) {
            break;
          } else {
            int al_index_name_symbol;
            al_index_name_symbol = dsat_acc2828dinitial_mpg78062[mappings_full_idx_symbol];
            if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
              break;
            } else {
              dsat_acc2828dinitial_der78106[((24 * ((j78189 * 3) + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dsat_acc2828dinitial_inv_mpg78063[((24 * ((j78189 * 3) + slice_idx)) + al_index_name_symbol)] >= 0) ? dsat_acc2828dinitial_der78106[((24 * ((j78189 * 3) + slice_idx)) + dsat_acc2828dinitial_inv_mpg78063[((24 * ((j78189 * 3) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) - (((((((((al_index_name_symbol < 24) ? ((dr4219dinitial_inv_mpg78067[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dr4219dinitial_der78187[((24 * (0 + slice_idx)) + dr4219dinitial_inv_mpg78067[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * sat_gms78046[i78185]) * dist278120) - (((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0) * (r78188[(0 + slice_idx)] * sat_gms78046[i78185]))) / (dist278120 * dist278120)) * sqrtl(dist278120)) - ((0.5 * (sqrtl((1.0 / dist278120)) * ((((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1) >= 0) ? ddist23094dinitial_der78119[((al_index_name_symbol < 24) ? ddist23094dinitial_inv_mpg78095[al_index_name_symbol] : -1)] : 0.0))) * ((r78188[(0 + slice_idx)] * sat_gms78046[i78185]) / dist278120))) / (sqrtl(dist278120) * sqrtl(dist278120))));
            }
          }
        }
      }
      for (int slice_idx = 0; slice_idx < (((j78189 * 3) + 3) - (j78189 * 3)); slice_idx++) {
        sat_acc78107[((j78189 * 3) + slice_idx)] = (sat_acc78107[((j78189 * 3) + slice_idx)] - (((sat_gms78046[i78185] * r78188[(0 + slice_idx)]) / dist278120) / sqrtl(dist278120)));
      }
    }
  }
  long double grav_acc_local78209[3] = { 0.0 };
  long double dgrav_acc_local4544dinitial_der78208[72] = { 0.0 };
  long double grav_acc78211[3] = { 0.0 };
  long double dgrav_acc4569dinitial_der78210[72] = { 0.0 };
  long double acc78213[3] = { 0.0 };
  long double dacc4588dinitial_der78212[72] = { 0.0 };
  long double rot78214[9] = { 0.0 };
  long double v_078216[7] = { 0.0 };
  long double dv_04620dinitial_der78215[168] = { 0.0 };
  long double dv_0_dx78218[7] = { 0.0 };
  long double ddv_0_dx4638dinitial_der78217[168] = { 0.0 };
  long double dv_0_dy78220[7] = { 0.0 };
  long double ddv_0_dy4660dinitial_der78219[168] = { 0.0 };
  long double dv_0_dz78222[7] = { 0.0 };
  long double ddv_0_dz4682dinitial_der78221[168] = { 0.0 };
  {
    long double jupiter_rotation_matrix78231[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix78231, t78041);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot78214[(0 + slice_idx)] = jupiter_rotation_matrix78231[slice_idx];
    }
  }
  for (int i78232 = 0; i78232 < 4; i78232++) {
    long double x78235;
    long double dx4794dinitial_der78234[24] = { 0.0 };
    for (int loop_var78239 = 0; loop_var78239 < 24; loop_var78239++) {
      int mapped_idx;
      mapped_idx = loop_var78239;
      int al_index_name_symbol;
      al_index_name_symbol = dx4794dinitial_mpg78082[loop_var78239];
      dx4794dinitial_der78234[mapped_idx] = (((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 0)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((i78232 * 6) + 0)) + dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 0)) + al_index_name_symbol)])] : 0.0) : 0.0) / 0.0004778945025452157572e0);
    }
    x78235 = (state78112[((i78232 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y78242;
    long double dy4843dinitial_der78241[24] = { 0.0 };
    for (int loop_var78246 = 0; loop_var78246 < 24; loop_var78246++) {
      int mapped_idx;
      mapped_idx = loop_var78246;
      int al_index_name_symbol;
      al_index_name_symbol = dy4843dinitial_mpg78098[loop_var78246];
      dy4843dinitial_der78241[mapped_idx] = (((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 1)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((i78232 * 6) + 1)) + dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 1)) + al_index_name_symbol)])] : 0.0) : 0.0) / 0.0004778945025452157572e0);
    }
    y78242 = (state78112[((i78232 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z78249;
    long double dz4892dinitial_der78248[24] = { 0.0 };
    for (int loop_var78253 = 0; loop_var78253 < 24; loop_var78253++) {
      int mapped_idx;
      mapped_idx = loop_var78253;
      int al_index_name_symbol;
      al_index_name_symbol = dz4892dinitial_mpg78096[loop_var78253];
      dz4892dinitial_der78248[mapped_idx] = (((al_index_name_symbol < 24) ? ((dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 2)) + al_index_name_symbol)] >= 0) ? dstate2924dinitial_der78111[((24 * ((i78232 * 6) + 2)) + dstate2924dinitial_inv_mpg78053[((24 * ((i78232 * 6) + 2)) + al_index_name_symbol)])] : 0.0) : 0.0) / 0.0004778945025452157572e0);
    }
    z78249 = (state78112[((i78232 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x78256;
    long double d_x4942dinitial_der78255[24] = { 0.0 };
    for (int loop_var78260 = 0; loop_var78260 < 24; loop_var78260++) {
      int mapped_idx;
      mapped_idx = loop_var78260;
      int al_index_name_symbol;
      al_index_name_symbol = d_x4942dinitial_mpg78076[loop_var78260];
      d_x4942dinitial_der78255[mapped_idx] = (((((((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1) >= 0) ? dx4794dinitial_der78234[((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1)] : 0.0) * rot78214[0]) + (((((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1) >= 0) ? dy4843dinitial_der78241[((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1)] : 0.0) * rot78214[3])) + (((((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1) >= 0) ? dz4892dinitial_der78248[((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1)] : 0.0) * rot78214[6]));
    }
    _x78256 = (((rot78214[0] * x78235) + (rot78214[3] * y78242)) + (rot78214[6] * z78249));
    long double _y78263;
    long double d_y4993dinitial_der78262[24] = { 0.0 };
    for (int loop_var78267 = 0; loop_var78267 < 24; loop_var78267++) {
      int mapped_idx;
      mapped_idx = loop_var78267;
      int al_index_name_symbol;
      al_index_name_symbol = d_y4993dinitial_mpg78102[loop_var78267];
      d_y4993dinitial_der78262[mapped_idx] = (((((((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1) >= 0) ? dx4794dinitial_der78234[((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1)] : 0.0) * rot78214[1]) + (((((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1) >= 0) ? dy4843dinitial_der78241[((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1)] : 0.0) * rot78214[4])) + (((((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1) >= 0) ? dz4892dinitial_der78248[((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1)] : 0.0) * rot78214[7]));
    }
    _y78263 = (((rot78214[1] * x78235) + (rot78214[4] * y78242)) + (rot78214[7] * z78249));
    long double _z78270;
    long double d_z5044dinitial_der78269[24] = { 0.0 };
    for (int loop_var78274 = 0; loop_var78274 < 24; loop_var78274++) {
      int mapped_idx;
      mapped_idx = loop_var78274;
      int al_index_name_symbol;
      al_index_name_symbol = d_z5044dinitial_mpg78056[loop_var78274];
      d_z5044dinitial_der78269[mapped_idx] = (((((((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1) >= 0) ? dx4794dinitial_der78234[((al_index_name_symbol < 24) ? dx4794dinitial_inv_mpg78083[al_index_name_symbol] : -1)] : 0.0) * rot78214[2]) + (((((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1) >= 0) ? dy4843dinitial_der78241[((al_index_name_symbol < 24) ? dy4843dinitial_inv_mpg78099[al_index_name_symbol] : -1)] : 0.0) * rot78214[5])) + (((((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1) >= 0) ? dz4892dinitial_der78248[((al_index_name_symbol < 24) ? dz4892dinitial_inv_mpg78097[al_index_name_symbol] : -1)] : 0.0) * rot78214[8]));
    }
    _z78270 = (((rot78214[2] * x78235) + (rot78214[5] * y78242)) + (rot78214[8] * z78249));
    long double _r278277;
    long double d_r25100dinitial_der78276[24] = { 0.0 };
    for (int loop_var78281 = 0; loop_var78281 < 24; loop_var78281++) {
      int mapped_idx;
      mapped_idx = loop_var78281;
      int al_index_name_symbol;
      al_index_name_symbol = d_r25100dinitial_mpg78104[loop_var78281];
      d_r25100dinitial_der78276[mapped_idx] = ((((_x78256 * ((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0)) + (_x78256 * ((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0))) + ((_y78263 * ((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0)) + (_y78263 * ((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0)))) + ((_z78270 * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0))));
    }
    _r278277 = (((_x78256 * _x78256) + (_y78263 * _y78263)) + (_z78270 * _z78270));
    long double _r78284;
    long double d_r5143dinitial_der78283[24] = { 0.0 };
    for (int loop_var78288 = 0; loop_var78288 < 24; loop_var78288++) {
      int mapped_idx;
      mapped_idx = loop_var78288;
      int al_index_name_symbol;
      al_index_name_symbol = d_r5143dinitial_mpg78060[loop_var78288];
      d_r5143dinitial_der78283[mapped_idx] = (0.5 * (sqrtl((1.0 / _r278277)) * ((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0)));
    }
    _r78284 = sqrtl(_r278277);
    long double _r378291;
    long double d_r35167dinitial_der78290[24] = { 0.0 };
    for (int loop_var78295 = 0; loop_var78295 < 24; loop_var78295++) {
      int mapped_idx;
      mapped_idx = loop_var78295;
      int al_index_name_symbol;
      al_index_name_symbol = d_r35167dinitial_mpg78070[loop_var78295];
      d_r35167dinitial_der78290[mapped_idx] = ((_r78284 * ((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0)) + (_r278277 * ((((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1) >= 0) ? d_r5143dinitial_der78283[((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1)] : 0.0)));
    }
    _r378291 = (_r78284 * _r278277);
    long double _r478298;
    long double d_r45191dinitial_der78297[24] = { 0.0 };
    for (int loop_var78302 = 0; loop_var78302 < 24; loop_var78302++) {
      int mapped_idx;
      mapped_idx = loop_var78302;
      int al_index_name_symbol;
      al_index_name_symbol = d_r45191dinitial_mpg78092[loop_var78302];
      d_r45191dinitial_der78297[mapped_idx] = ((_r278277 * ((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0)) + (_r278277 * ((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0)));
    }
    _r478298 = (_r278277 * _r278277);
    long double _r578305;
    long double d_r55216dinitial_der78304[24] = { 0.0 };
    for (int loop_var78309 = 0; loop_var78309 < 24; loop_var78309++) {
      int mapped_idx;
      mapped_idx = loop_var78309;
      int al_index_name_symbol;
      al_index_name_symbol = d_r55216dinitial_mpg78090[loop_var78309];
      d_r55216dinitial_der78304[mapped_idx] = ((_r478298 * ((((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1) >= 0) ? d_r5143dinitial_der78283[((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1)] : 0.0)) + (_r78284 * ((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0)));
    }
    _r578305 = (_r478298 * _r78284);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg78072[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der78215[((24 * 0) + mapped_idx)] = (-1.0 * (1.0e0 * ((1.0 / (_r78284 * _r78284)) * ((((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1) >= 0) ? d_r5143dinitial_der78283[((al_index_name_symbol < 24) ? d_r5143dinitial_inv_mpg78061[al_index_name_symbol] : -1)] : 0.0))));
        }
      }
    }
    v_078216[0] = (1.0e0 / _r78284);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dv_04620dinitial_mpg78072[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dv_04620dinitial_der78215[((24 * 1) + mapped_idx)] = ((((((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * 1.7320508075688772936e0) * _r378291) - (((((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1) >= 0) ? d_r35167dinitial_der78290[((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1)] : 0.0) * (_z78270 * 1.7320508075688772936e0))) / (_r378291 * _r378291));
        }
      }
    }
    v_078216[1] = ((1.7320508075688772936e0 * _z78270) / _r378291);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg78058[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der78217[((24 * 0) + mapped_idx)] = (((-((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0) * _r378291) - (((((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1) >= 0) ? d_r35167dinitial_der78290[((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1)] : 0.0) * -_x78256)) / (_r378291 * _r378291));
        }
      }
    }
    dv_0_dx78218[0] = (-_x78256 / _r378291);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg78080[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der78219[((24 * 0) + mapped_idx)] = (((-((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0) * _r378291) - (((((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1) >= 0) ? d_r35167dinitial_der78290[((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1)] : 0.0) * -_y78263)) / (_r378291 * _r378291));
        }
      }
    }
    dv_0_dy78220[0] = (-_y78263 / _r378291);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg78064[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der78221[((24 * 0) + mapped_idx)] = (((-((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * _r378291) - (((((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1) >= 0) ? d_r35167dinitial_der78290[((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1)] : 0.0) * -_z78270)) / (_r378291 * _r378291));
        }
      }
    }
    dv_0_dz78222[0] = (-_z78270 / _r378291);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dx4638dinitial_mpg78058[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dx4638dinitial_der78217[((24 * 1) + mapped_idx)] = ((((((_z78270 * -5.1961524227066318805e0) * ((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0)) + (_x78256 * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * -5.1961524227066318805e0))) * _r578305) - (((((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1) >= 0) ? d_r55216dinitial_der78304[((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1)] : 0.0) * ((_z78270 * -5.1961524227066318805e0) * _x78256))) / (_r578305 * _r578305));
        }
      }
    }
    dv_0_dx78218[1] = (((-5.1961524227066318805e0 * _z78270) * _x78256) / _r578305);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dy4660dinitial_mpg78080[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dy4660dinitial_der78219[((24 * 1) + mapped_idx)] = ((((((_z78270 * -5.1961524227066318805e0) * ((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0)) + (_y78263 * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * -5.1961524227066318805e0))) * _r578305) - (((((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1) >= 0) ? d_r55216dinitial_der78304[((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1)] : 0.0) * ((_z78270 * -5.1961524227066318805e0) * _y78263))) / (_r578305 * _r578305));
        }
      }
    }
    dv_0_dy78220[1] = (((-5.1961524227066318805e0 * _z78270) * _y78263) / _r578305);
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 168)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = ddv_0_dz4682dinitial_mpg78064[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          ddv_0_dz4682dinitial_der78221[((24 * 1) + mapped_idx)] = (((-1.0 * (1.0e0 * ((1.0 / (_r378291 * _r378291)) * ((((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1) >= 0) ? d_r35167dinitial_der78290[((al_index_name_symbol < 24) ? d_r35167dinitial_inv_mpg78071[al_index_name_symbol] : -1)] : 0.0)))) - ((((((_z78270 * 3.0e0) * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * 3.0e0))) * _r578305) - (((((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1) >= 0) ? d_r55216dinitial_der78304[((al_index_name_symbol < 24) ? d_r55216dinitial_inv_mpg78091[al_index_name_symbol] : -1)] : 0.0) * ((_z78270 * 3.0e0) * _z78270))) / (_r578305 * _r578305))) * 1.7320508075688772936e0);
        }
      }
    }
    dv_0_dz78222[1] = (1.7320508075688772936e0 * ((1.0e0 / _r378291) - (((3.0e0 * _z78270) * _z78270) / _r578305)));
    for (int n78343 = 2; n78343 < 7; n78343++) {
      long double coef178345;
      coef178345 = sqrtl(((((2.0e0 * ((long double) n78343)) - 1.0e0) * ((long double) ((2 * n78343) + 1))) / ((long double) (n78343 * n78343))));
      long double coef278349;
      coef278349 = sqrtl(((((((long double) n78343) - 1.0e0) * ((long double) (n78343 - 1))) * ((long double) ((2 * n78343) + 1))) / ((long double) ((n78343 * n78343) * ((2 * n78343) - 3)))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n78343));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dv_04620dinitial_mpg78072[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dv_04620dinitial_der78215[((24 * n78343) + mapped_idx)] = (((((((v_078216[(n78343 - 1)] * coef178345) * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 1)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0) * coef178345))) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 1)] * coef178345) * _z78270))) / (_r278277 * _r278277)) - ((((((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 2)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * coef278349) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * (v_078216[(n78343 - 2)] * coef278349))) / (_r278277 * _r278277)));
          }
        }
      }
      v_078216[n78343] = ((((coef178345 * v_078216[(n78343 - 1)]) * _z78270) / _r278277) - ((coef278349 * v_078216[(n78343 - 2)]) / _r278277));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n78343));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dx4638dinitial_mpg78058[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dx4638dinitial_der78217[((24 * n78343) + mapped_idx)] = ((((_z78270 * coef178345) * ((((((al_index_name_symbol < 24) ? ((ddv_0_dx4638dinitial_inv_mpg78059[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? ddv_0_dx4638dinitial_der78217[((24 * (n78343 - 1)) + ddv_0_dx4638dinitial_inv_mpg78059[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * dv_0_dx78218[(n78343 - 1)])) / (_r278277 * _r278277)) - ((((((v_078216[(n78343 - 1)] * 2.0e0) * ((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0)) + (_x78256 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 1)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 1)] * 2.0e0) * _x78256))) / (_r478298 * _r478298)))) + (((dv_0_dx78218[(n78343 - 1)] / _r278277) - (((v_078216[(n78343 - 1)] * 2.0e0) * _x78256) / _r478298)) * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * coef178345))) - (((((((al_index_name_symbol < 24) ? ((ddv_0_dx4638dinitial_inv_mpg78059[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? ddv_0_dx4638dinitial_der78217[((24 * (n78343 - 2)) + ddv_0_dx4638dinitial_inv_mpg78059[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * dv_0_dx78218[(n78343 - 2)])) / (_r278277 * _r278277)) - ((((((v_078216[(n78343 - 2)] * 2.0e0) * ((((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1) >= 0) ? d_x4942dinitial_der78255[((al_index_name_symbol < 24) ? d_x4942dinitial_inv_mpg78077[al_index_name_symbol] : -1)] : 0.0)) + (_x78256 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 2)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 2)] * 2.0e0) * _x78256))) / (_r478298 * _r478298))) * coef278349));
          }
        }
      }
      dv_0_dx78218[n78343] = (((coef178345 * _z78270) * ((dv_0_dx78218[(n78343 - 1)] / _r278277) - (((v_078216[(n78343 - 1)] * 2.0e0) * _x78256) / _r478298))) - (coef278349 * ((dv_0_dx78218[(n78343 - 2)] / _r278277) - (((v_078216[(n78343 - 2)] * 2.0e0) * _x78256) / _r478298))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n78343));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dy4660dinitial_mpg78080[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dy4660dinitial_der78219[((24 * n78343) + mapped_idx)] = ((((_z78270 * coef178345) * ((((((al_index_name_symbol < 24) ? ((ddv_0_dy4660dinitial_inv_mpg78081[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? ddv_0_dy4660dinitial_der78219[((24 * (n78343 - 1)) + ddv_0_dy4660dinitial_inv_mpg78081[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * dv_0_dy78220[(n78343 - 1)])) / (_r278277 * _r278277)) - ((((((v_078216[(n78343 - 1)] * 2.0e0) * ((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0)) + (_y78263 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 1)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 1)] * 2.0e0) * _y78263))) / (_r478298 * _r478298)))) + (((dv_0_dy78220[(n78343 - 1)] / _r278277) - (((v_078216[(n78343 - 1)] * 2.0e0) * _y78263) / _r478298)) * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * coef178345))) - (((((((al_index_name_symbol < 24) ? ((ddv_0_dy4660dinitial_inv_mpg78081[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? ddv_0_dy4660dinitial_der78219[((24 * (n78343 - 2)) + ddv_0_dy4660dinitial_inv_mpg78081[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * dv_0_dy78220[(n78343 - 2)])) / (_r278277 * _r278277)) - ((((((v_078216[(n78343 - 2)] * 2.0e0) * ((((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1) >= 0) ? d_y4993dinitial_der78262[((al_index_name_symbol < 24) ? d_y4993dinitial_inv_mpg78103[al_index_name_symbol] : -1)] : 0.0)) + (_y78263 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 2)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 2)] * 2.0e0) * _y78263))) / (_r478298 * _r478298))) * coef278349));
          }
        }
      }
      dv_0_dy78220[n78343] = (((coef178345 * _z78270) * ((dv_0_dy78220[(n78343 - 1)] / _r278277) - (((v_078216[(n78343 - 1)] * 2.0e0) * _y78263) / _r478298))) - (coef278349 * ((dv_0_dy78220[(n78343 - 2)] / _r278277) - (((v_078216[(n78343 - 2)] * 2.0e0) * _y78263) / _r478298))));
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * n78343));
        if ((mappings_full_idx_symbol >= 168)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = ddv_0_dz4682dinitial_mpg78064[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            ddv_0_dz4682dinitial_der78221[((24 * n78343) + mapped_idx)] = ((((((((dv_0_dz78222[(n78343 - 1)] * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * ((al_index_name_symbol < 24) ? ((ddv_0_dz4682dinitial_inv_mpg78065[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? ddv_0_dz4682dinitial_der78221[((24 * (n78343 - 1)) + ddv_0_dz4682dinitial_inv_mpg78065[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0))) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * (dv_0_dz78222[(n78343 - 1)] * _z78270))) / (_r278277 * _r278277)) + ((v_078216[(n78343 - 1)] * ((-1.0 * (1.0e0 * ((1.0 / (_r278277 * _r278277)) * ((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0)))) - ((((((_z78270 * 2.0e0) * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * (((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((_z78270 * 2.0e0) * _z78270))) / (_r478298 * _r478298)))) + (((1.0e0 / _r278277) - (((_z78270 * 2.0e0) * _z78270) / _r478298)) * ((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 1)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 1)) + al_index_name_symbol)])] : 0.0) : 0.0)))) * coef178345) - (((((((al_index_name_symbol < 24) ? ((ddv_0_dz4682dinitial_inv_mpg78065[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? ddv_0_dz4682dinitial_der78221[((24 * (n78343 - 2)) + ddv_0_dz4682dinitial_inv_mpg78065[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * _r278277) - (((((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1) >= 0) ? d_r25100dinitial_der78276[((al_index_name_symbol < 24) ? d_r25100dinitial_inv_mpg78105[al_index_name_symbol] : -1)] : 0.0) * dv_0_dz78222[(n78343 - 2)])) / (_r278277 * _r278277)) - ((((((v_078216[(n78343 - 2)] * 2.0e0) * ((((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1) >= 0) ? d_z5044dinitial_der78269[((al_index_name_symbol < 24) ? d_z5044dinitial_inv_mpg78057[al_index_name_symbol] : -1)] : 0.0)) + (_z78270 * (((al_index_name_symbol < 24) ? ((dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)] >= 0) ? dv_04620dinitial_der78215[((24 * (n78343 - 2)) + dv_04620dinitial_inv_mpg78073[((24 * (n78343 - 2)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.0e0))) * _r478298) - (((((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1) >= 0) ? d_r45191dinitial_der78297[((al_index_name_symbol < 24) ? d_r45191dinitial_inv_mpg78093[al_index_name_symbol] : -1)] : 0.0) * ((v_078216[(n78343 - 2)] * 2.0e0) * _z78270))) / (_r478298 * _r478298))) * coef278349));
          }
        }
      }
      dv_0_dz78222[n78343] = ((coef178345 * (((dv_0_dz78222[(n78343 - 1)] * _z78270) / _r278277) + (v_078216[(n78343 - 1)] * ((1.0e0 / _r278277) - (((2.0e0 * _z78270) * _z78270) / _r478298))))) - (coef278349 * ((dv_0_dz78222[(n78343 - 2)] / _r278277) - (((v_078216[(n78343 - 2)] * 2.0e0) * _z78270) / _r478298))));
    }
    long double dpotential_dx78370;
    long double ddpotential_dx6387dinitial_der78369[24] = { 0.0 };
    for (int loop_var78374 = 0; loop_var78374 < 24; loop_var78374++) {
      int mapped_idx;
      mapped_idx = loop_var78374;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dx6387dinitial_mpg78050[loop_var78374];
      ddpotential_dx6387dinitial_der78369[mapped_idx] = 0.0;
    }
    dpotential_dx78370 = ((long double) 0);
    long double dpotential_dy78376;
    long double ddpotential_dy6414dinitial_der78375[24] = { 0.0 };
    for (int loop_var78380 = 0; loop_var78380 < 24; loop_var78380++) {
      int mapped_idx;
      mapped_idx = loop_var78380;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dy6414dinitial_mpg78088[loop_var78380];
      ddpotential_dy6414dinitial_der78375[mapped_idx] = 0.0;
    }
    dpotential_dy78376 = ((long double) 0);
    long double dpotential_dz78382;
    long double ddpotential_dz6441dinitial_der78381[24] = { 0.0 };
    for (int loop_var78386 = 0; loop_var78386 < 24; loop_var78386++) {
      int mapped_idx;
      mapped_idx = loop_var78386;
      int al_index_name_symbol;
      al_index_name_symbol = ddpotential_dz6441dinitial_mpg78068[loop_var78386];
      ddpotential_dz6441dinitial_der78381[mapped_idx] = 0.0;
    }
    dpotential_dz78382 = ((long double) 0);
    for (int n78387 = 2; n78387 < 7; n78387++) {
      for (int loop_var78392 = 0; loop_var78392 < 24; loop_var78392++) {
        int mapped_idx;
        mapped_idx = loop_var78392;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dx6387dinitial_mpg78050[loop_var78392];
        ddpotential_dx6387dinitial_der78369[mapped_idx] = (((((al_index_name_symbol < 24) ? ddpotential_dx6387dinitial_inv_mpg78051[al_index_name_symbol] : -1) >= 0) ? ddpotential_dx6387dinitial_der78369[((al_index_name_symbol < 24) ? ddpotential_dx6387dinitial_inv_mpg78051[al_index_name_symbol] : -1)] : 0.0) + (((al_index_name_symbol < 24) ? ((ddv_0_dx4638dinitial_inv_mpg78059[((24 * n78387) + al_index_name_symbol)] >= 0) ? ddv_0_dx4638dinitial_der78217[((24 * n78387) + ddv_0_dx4638dinitial_inv_mpg78059[((24 * n78387) + al_index_name_symbol)])] : 0.0) : 0.0) * central_grav77975[n78387]));
      }
      dpotential_dx78370 = (dpotential_dx78370 + (dv_0_dx78218[n78387] * central_grav77975[n78387]));
      for (int loop_var78397 = 0; loop_var78397 < 24; loop_var78397++) {
        int mapped_idx;
        mapped_idx = loop_var78397;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dy6414dinitial_mpg78088[loop_var78397];
        ddpotential_dy6414dinitial_der78375[mapped_idx] = (((((al_index_name_symbol < 24) ? ddpotential_dy6414dinitial_inv_mpg78089[al_index_name_symbol] : -1) >= 0) ? ddpotential_dy6414dinitial_der78375[((al_index_name_symbol < 24) ? ddpotential_dy6414dinitial_inv_mpg78089[al_index_name_symbol] : -1)] : 0.0) + (((al_index_name_symbol < 24) ? ((ddv_0_dy4660dinitial_inv_mpg78081[((24 * n78387) + al_index_name_symbol)] >= 0) ? ddv_0_dy4660dinitial_der78219[((24 * n78387) + ddv_0_dy4660dinitial_inv_mpg78081[((24 * n78387) + al_index_name_symbol)])] : 0.0) : 0.0) * central_grav77975[n78387]));
      }
      dpotential_dy78376 = (dpotential_dy78376 + (dv_0_dy78220[n78387] * central_grav77975[n78387]));
      for (int loop_var78402 = 0; loop_var78402 < 24; loop_var78402++) {
        int mapped_idx;
        mapped_idx = loop_var78402;
        int al_index_name_symbol;
        al_index_name_symbol = ddpotential_dz6441dinitial_mpg78068[loop_var78402];
        ddpotential_dz6441dinitial_der78381[mapped_idx] = (((((al_index_name_symbol < 24) ? ddpotential_dz6441dinitial_inv_mpg78069[al_index_name_symbol] : -1) >= 0) ? ddpotential_dz6441dinitial_der78381[((al_index_name_symbol < 24) ? ddpotential_dz6441dinitial_inv_mpg78069[al_index_name_symbol] : -1)] : 0.0) + (((al_index_name_symbol < 24) ? ((ddv_0_dz4682dinitial_inv_mpg78065[((24 * n78387) + al_index_name_symbol)] >= 0) ? ddv_0_dz4682dinitial_der78221[((24 * n78387) + ddv_0_dz4682dinitial_inv_mpg78065[((24 * n78387) + al_index_name_symbol)])] : 0.0) : 0.0) * central_grav77975[n78387]));
      }
      dpotential_dz78382 = (dpotential_dz78382 + (dv_0_dz78222[n78387] * central_grav77975[n78387]));
    }
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 0));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg78074[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der78208[((24 * 0) + mapped_idx)] = ((((al_index_name_symbol < 24) ? ddpotential_dx6387dinitial_inv_mpg78051[al_index_name_symbol] : -1) >= 0) ? ddpotential_dx6387dinitial_der78369[((al_index_name_symbol < 24) ? ddpotential_dx6387dinitial_inv_mpg78051[al_index_name_symbol] : -1)] : 0.0);
        }
      }
    }
    grav_acc_local78209[0] = dpotential_dx78370;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 1));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg78074[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der78208[((24 * 1) + mapped_idx)] = ((((al_index_name_symbol < 24) ? ddpotential_dy6414dinitial_inv_mpg78089[al_index_name_symbol] : -1) >= 0) ? ddpotential_dy6414dinitial_der78375[((al_index_name_symbol < 24) ? ddpotential_dy6414dinitial_inv_mpg78089[al_index_name_symbol] : -1)] : 0.0);
        }
      }
    }
    grav_acc_local78209[1] = dpotential_dy78376;
    for (int mapped_idx = 0;; mapped_idx++) {
      int mappings_full_idx_symbol;
      mappings_full_idx_symbol = (mapped_idx + (24 * 2));
      if ((mappings_full_idx_symbol >= 72)) {
        break;
      } else {
        int al_index_name_symbol;
        al_index_name_symbol = dgrav_acc_local4544dinitial_mpg78074[mappings_full_idx_symbol];
        if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
          break;
        } else {
          dgrav_acc_local4544dinitial_der78208[((24 * 2) + mapped_idx)] = ((((al_index_name_symbol < 24) ? ddpotential_dz6441dinitial_inv_mpg78069[al_index_name_symbol] : -1) >= 0) ? ddpotential_dz6441dinitial_der78381[((al_index_name_symbol < 24) ? ddpotential_dz6441dinitial_inv_mpg78069[al_index_name_symbol] : -1)] : 0.0);
        }
      }
    }
    grav_acc_local78209[2] = dpotential_dz78382;
    for (int k78416 = 0; k78416 < 3; k78416++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * k78416));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dgrav_acc4569dinitial_mpg78078[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dgrav_acc4569dinitial_der78210[((24 * k78416) + mapped_idx)] = ((((((al_index_name_symbol < 24) ? ((dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 0) + al_index_name_symbol)] >= 0) ? dgrav_acc_local4544dinitial_der78208[((24 * 0) + dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 0) + al_index_name_symbol)])] : 0.0) : 0.0) * rot78214[((3 * k78416) + 0)]) + (((al_index_name_symbol < 24) ? ((dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 1) + al_index_name_symbol)] >= 0) ? dgrav_acc_local4544dinitial_der78208[((24 * 1) + dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 1) + al_index_name_symbol)])] : 0.0) : 0.0) * rot78214[((3 * k78416) + 1)])) + (((al_index_name_symbol < 24) ? ((dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 2) + al_index_name_symbol)] >= 0) ? dgrav_acc_local4544dinitial_der78208[((24 * 2) + dgrav_acc_local4544dinitial_inv_mpg78075[((24 * 2) + al_index_name_symbol)])] : 0.0) : 0.0) * rot78214[((3 * k78416) + 2)])) / 2.2838315556293922983e-07);
          }
        }
      }
      grav_acc78211[k78416] = ((((rot78214[((3 * k78416) + 0)] * grav_acc_local78209[0]) + (rot78214[((3 * k78416) + 1)] * grav_acc_local78209[1])) + (rot78214[((3 * k78416) + 2)] * grav_acc_local78209[2])) / 2.2838315556293922983e-07);
    }
    for (int slice_idx = 0; slice_idx < (((i78232 * 3) + 3) - (i78232 * 3)); slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * ((i78232 * 3) + slice_idx)));
        if ((mappings_full_idx_symbol >= 288)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dsat_acc2828dinitial_mpg78062[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dsat_acc2828dinitial_der78106[((24 * ((i78232 * 3) + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78232 * 3) + slice_idx)) + al_index_name_symbol)] >= 0) ? dsat_acc2828dinitial_der78106[((24 * ((i78232 * 3) + slice_idx)) + dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78232 * 3) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) + (((al_index_name_symbol < 24) ? ((dgrav_acc4569dinitial_inv_mpg78079[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dgrav_acc4569dinitial_der78210[((24 * (0 + slice_idx)) + dgrav_acc4569dinitial_inv_mpg78079[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * 2.8247609439046209905e-07));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < (((i78232 * 3) + 3) - (i78232 * 3)); slice_idx++) {
      sat_acc78107[((i78232 * 3) + slice_idx)] = (sat_acc78107[((i78232 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc78211[(0 + slice_idx)]));
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dcentral_acc2846dinitial_mpg78054[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) - (((al_index_name_symbol < 24) ? ((dgrav_acc4569dinitial_inv_mpg78079[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dgrav_acc4569dinitial_der78210[((24 * (0 + slice_idx)) + dgrav_acc4569dinitial_inv_mpg78079[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) * sat_gms78046[i78232]));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      central_acc78109[(0 + slice_idx)] = (central_acc78109[(0 + slice_idx)] - (sat_gms78046[i78232] * grav_acc78211[(0 + slice_idx)]));
    }
  }
  for (int i78430 = 0; i78430 < 4; i78430++) {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      for (int mapped_idx = 0;; mapped_idx++) {
        int mappings_full_idx_symbol;
        mappings_full_idx_symbol = (mapped_idx + (24 * (0 + slice_idx)));
        if ((mappings_full_idx_symbol >= 72)) {
          break;
        } else {
          int al_index_name_symbol;
          al_index_name_symbol = dacc4588dinitial_mpg78100[mappings_full_idx_symbol];
          if (((al_index_name_symbol < 0) || (mapped_idx >= 24))) {
            break;
          } else {
            dacc4588dinitial_der78212[((24 * (0 + slice_idx)) + mapped_idx)] = (((al_index_name_symbol < 24) ? ((dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78430 * 3) + slice_idx)) + al_index_name_symbol)] >= 0) ? dsat_acc2828dinitial_der78106[((24 * ((i78430 * 3) + slice_idx)) + dsat_acc2828dinitial_inv_mpg78063[((24 * ((i78430 * 3) + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0) - ((al_index_name_symbol < 24) ? ((dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)] >= 0) ? dcentral_acc2846dinitial_der78108[((24 * (0 + slice_idx)) + dcentral_acc2846dinitial_inv_mpg78055[((24 * (0 + slice_idx)) + al_index_name_symbol)])] : 0.0) : 0.0));
          }
        }
      }
    }
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      acc78213[(0 + slice_idx)] = (sat_acc78107[((i78430 * 3) + slice_idx)] - central_acc78109[(0 + slice_idx)]);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i78430 * 6) + 3) - (i78430 * 6)); slice_idx++) {
        jupsatsystem78040[((i78430 * 6) + slice_idx)] = state78112[(((i78430 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i78430 * 6) + 6) - ((i78430 * 6) + 3)); slice_idx++) {
        jupsatsystem78040[(((i78430 * 6) + 3) + slice_idx)] = acc78213[(0 + slice_idx)];
      }
    }
    for (int j78442 = 0; j78442 < 3; j78442++) {
      {
        for (int slice_idx = 0; slice_idx < ((24 + (((((i78430 * 6) + j78442) + 1) * 4) * 6)) - (24 + ((((i78430 * 6) + j78442) * 4) * 6))); slice_idx++) {
          jupsatsystem78040[((24 + ((((i78430 * 6) + j78442) * 4) * 6)) + slice_idx)] = state_and_derivatives78042[((24 + (((((i78430 * 6) + j78442) + 3) * 4) * 6)) + slice_idx)];
        }
      }
      for (int slice_idx = 0; slice_idx < ((24 + (((((i78430 * 6) + j78442) + 4) * 4) * 6)) - (24 + (((((i78430 * 6) + j78442) + 3) * 4) * 6))); slice_idx++) {
        int al_index_name_symbol;
        al_index_name_symbol = (slice_idx + 0);
        int func_slice_idx;
        func_slice_idx = (slice_idx + (24 + (((((i78430 * 6) + j78442) + 3) * 4) * 6)));
        jupsatsystem78040[func_slice_idx] = ((al_index_name_symbol < 24) ? ((dacc4588dinitial_inv_mpg78101[((24 * j78442) + al_index_name_symbol)] >= 0) ? dacc4588dinitial_der78212[((24 * j78442) + dacc4588dinitial_inv_mpg78101[((24 * j78442) + al_index_name_symbol)])] : 0.0) : 0.0);
      }
    }
  }
  return 0;
}

int jupsatsystem_noderiv(long double *restrict jupsatsystem_noderiv78449, long double t78450, long double *restrict state78451, long double *restrict central_pos78452, long double *restrict perturb_gms78453, long double *restrict perturb_pos78454, long double *restrict sat_gms78455) {
  long double sat_acc78457[12] = { 0.0 };
  long double central_acc78458[3] = { 0.0 };
  long double dist278459;
  long double dist378460;
  for (int i78461 = 0; i78461 < 4; i78461++) {
    long double r78463[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r78463[(0 + slice_idx)] = state78451[((i78461 * 6) + slice_idx)];
      }
    }
    dist278459 = (((r78463[0] * r78463[0]) + (r78463[1] * r78463[1])) + (r78463[2] * r78463[2]));
    dist378460 = (dist278459 * sqrtl(dist278459));
    {
      for (int slice_idx = 0; slice_idx < (((i78461 * 3) + 3) - (i78461 * 3)); slice_idx++) {
        sat_acc78457[((i78461 * 3) + slice_idx)] = ((-2.8247609439046209905e-07 * r78463[(0 + slice_idx)]) / dist378460);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc78458[(0 + slice_idx)] = (central_acc78458[(0 + slice_idx)] + ((sat_gms78455[i78461] * r78463[(0 + slice_idx)]) / dist378460));
      }
    }
  }
  for (int i78479 = 0; i78479 < 4; i78479++) {
    long double r78481[3] = { 0.0 };
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        r78481[(0 + slice_idx)] = (perturb_pos78454[((i78479 * 3) + slice_idx)] - central_pos78452[(0 + slice_idx)]);
      }
    }
    dist278459 = (((r78481[0] * r78481[0]) + (r78481[1] * r78481[1])) + (r78481[2] * r78481[2]));
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc78458[(0 + slice_idx)] = (central_acc78458[(0 + slice_idx)] + (((perturb_gms78453[i78479] * r78481[(0 + slice_idx)]) / dist278459) / sqrtl(dist278459)));
      }
    }
    for (int j78491 = 0; j78491 < 4; j78491++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r78481[(0 + slice_idx)] = (perturb_pos78454[((i78479 * 3) + slice_idx)] - (state78451[((j78491 * 6) + slice_idx)] + central_pos78452[(0 + slice_idx)]));
        }
      }
      dist278459 = (((r78481[0] * r78481[0]) + (r78481[1] * r78481[1])) + (r78481[2] * r78481[2]));
      {
        for (int slice_idx = 0; slice_idx < (((j78491 * 3) + 3) - (j78491 * 3)); slice_idx++) {
          sat_acc78457[((j78491 * 3) + slice_idx)] = (sat_acc78457[((j78491 * 3) + slice_idx)] + (((perturb_gms78453[i78479] * r78481[(0 + slice_idx)]) / dist278459) / sqrtl(dist278459)));
        }
      }
    }
  }
  for (int i78502 = 1; i78502 < 4; i78502++) {
    long double r78504[3] = { 0.0 };
    for (int j78505 = 0; j78505 < i78502; j78505++) {
      {
        for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
          r78504[(0 + slice_idx)] = (state78451[((j78505 * 6) + slice_idx)] - state78451[((i78502 * 6) + slice_idx)]);
        }
      }
      dist278459 = (((r78504[0] * r78504[0]) + (r78504[1] * r78504[1])) + (r78504[2] * r78504[2]));
      {
        for (int slice_idx = 0; slice_idx < (((i78502 * 3) + 3) - (i78502 * 3)); slice_idx++) {
          sat_acc78457[((i78502 * 3) + slice_idx)] = (sat_acc78457[((i78502 * 3) + slice_idx)] + (((sat_gms78455[j78505] * r78504[(0 + slice_idx)]) / dist278459) / sqrtl(dist278459)));
        }
      }
      {
        for (int slice_idx = 0; slice_idx < (((j78505 * 3) + 3) - (j78505 * 3)); slice_idx++) {
          sat_acc78457[((j78505 * 3) + slice_idx)] = (sat_acc78457[((j78505 * 3) + slice_idx)] - (((sat_gms78455[i78502] * r78504[(0 + slice_idx)]) / dist278459) / sqrtl(dist278459)));
        }
      }
    }
  }
  long double grav_acc_local78519[3] = { 0.0 };
  long double grav_acc78520[3] = { 0.0 };
  long double acc78521[3] = { 0.0 };
  long double rot78522[9] = { 0.0 };
  long double v_078523[7] = { 0.0 };
  long double dv_0_dx78524[7] = { 0.0 };
  long double dv_0_dy78525[7] = { 0.0 };
  long double dv_0_dz78526[7] = { 0.0 };
  {
    long double jupiter_rotation_matrix78535[9] = { 0.0 };
    jupiter_rotation_matrix(jupiter_rotation_matrix78535, t78450);
    for (int slice_idx = 0; slice_idx < 9; slice_idx++) {
      rot78522[(0 + slice_idx)] = jupiter_rotation_matrix78535[slice_idx];
    }
  }
  for (int i78536 = 0; i78536 < 4; i78536++) {
    long double x78538;
    x78538 = (state78451[((i78536 * 6) + 0)] / 0.0004778945025452157572e0);
    long double y78542;
    y78542 = (state78451[((i78536 * 6) + 1)] / 0.0004778945025452157572e0);
    long double z78546;
    z78546 = (state78451[((i78536 * 6) + 2)] / 0.0004778945025452157572e0);
    long double _x78550;
    _x78550 = (((rot78522[0] * x78538) + (rot78522[3] * y78542)) + (rot78522[6] * z78546));
    long double _y78554;
    _y78554 = (((rot78522[1] * x78538) + (rot78522[4] * y78542)) + (rot78522[7] * z78546));
    long double _z78558;
    _z78558 = (((rot78522[2] * x78538) + (rot78522[5] * y78542)) + (rot78522[8] * z78546));
    long double _r278562;
    _r278562 = (((_x78550 * _x78550) + (_y78554 * _y78554)) + (_z78558 * _z78558));
    long double r78566;
    r78566 = sqrtl(_r278562);
    long double _r378570;
    _r378570 = (r78566 * _r278562);
    long double _r478574;
    _r478574 = (_r278562 * _r278562);
    long double _r578578;
    _r578578 = (_r478574 * r78566);
    v_078523[0] = (1.0e0 / r78566);
    v_078523[1] = ((1.7320508075688772936e0 * _z78558) / _r378570);
    dv_0_dx78524[0] = (-_x78550 / _r378570);
    dv_0_dy78525[0] = (-_y78554 / _r378570);
    dv_0_dz78526[0] = (-_z78558 / _r378570);
    dv_0_dx78524[1] = (((-5.1961524227066318805e0 * _z78558) * _x78550) / _r578578);
    dv_0_dy78525[1] = (((-5.1961524227066318805e0 * _z78558) * _y78554) / _r578578);
    dv_0_dz78526[1] = (1.7320508075688772936e0 * ((1.0e0 / _r378570) - (((3.0e0 * _z78558) * _z78558) / _r578578)));
    for (int n78606 = 2; n78606 < 7; n78606++) {
      long double coef178608;
      coef178608 = sqrtl(((((2.0e0 * ((long double) n78606)) - 1.0e0) * ((long double) ((2 * n78606) + 1))) / ((long double) (n78606 * n78606))));
      long double coef278612;
      coef278612 = sqrtl(((((((long double) n78606) - 1.0e0) * ((long double) (n78606 - 1))) * ((long double) ((2 * n78606) + 1))) / ((long double) ((n78606 * n78606) * ((2 * n78606) - 3)))));
      v_078523[n78606] = ((((coef178608 * v_078523[(n78606 - 1)]) * _z78558) / _r278562) - ((coef278612 * v_078523[(n78606 - 2)]) / _r278562));
      dv_0_dx78524[n78606] = (((coef178608 * _z78558) * ((dv_0_dx78524[(n78606 - 1)] / _r278562) - (((v_078523[(n78606 - 1)] * 2.0e0) * _x78550) / _r478574))) - (coef278612 * ((dv_0_dx78524[(n78606 - 2)] / _r278562) - (((v_078523[(n78606 - 2)] * 2.0e0) * _x78550) / _r478574))));
      dv_0_dy78525[n78606] = (((coef178608 * _z78558) * ((dv_0_dy78525[(n78606 - 1)] / _r278562) - (((v_078523[(n78606 - 1)] * 2.0e0) * _y78554) / _r478574))) - (coef278612 * ((dv_0_dy78525[(n78606 - 2)] / _r278562) - (((v_078523[(n78606 - 2)] * 2.0e0) * _y78554) / _r478574))));
      dv_0_dz78526[n78606] = ((coef178608 * (((dv_0_dz78526[(n78606 - 1)] * _z78558) / _r278562) + (v_078523[(n78606 - 1)] * ((1.0e0 / _r278562) - (((2.0e0 * _z78558) * _z78558) / _r478574))))) - (coef278612 * ((dv_0_dz78526[(n78606 - 2)] / _r278562) - (((v_078523[(n78606 - 2)] * 2.0e0) * _z78558) / _r478574))));
    }
    long double dpotential_dx78628;
    dpotential_dx78628 = ((long double) 0);
    long double dpotential_dy78632;
    dpotential_dy78632 = ((long double) 0);
    long double dpotential_dz78636;
    dpotential_dz78636 = ((long double) 0);
    for (int n78640 = 2; n78640 < 7; n78640++) {
      dpotential_dx78628 = (dpotential_dx78628 + (dv_0_dx78524[n78640] * central_grav77975[n78640]));
      dpotential_dy78632 = (dpotential_dy78632 + (dv_0_dy78525[n78640] * central_grav77975[n78640]));
      dpotential_dz78636 = (dpotential_dz78636 + (dv_0_dz78526[n78640] * central_grav77975[n78640]));
    }
    grav_acc_local78519[0] = dpotential_dx78628;
    grav_acc_local78519[1] = dpotential_dy78632;
    grav_acc_local78519[2] = dpotential_dz78636;
    for (int k78660 = 0; k78660 < 3; k78660++) {
      grav_acc78520[k78660] = ((((rot78522[((3 * k78660) + 0)] * grav_acc_local78519[0]) + (rot78522[((3 * k78660) + 1)] * grav_acc_local78519[1])) + (rot78522[((3 * k78660) + 2)] * grav_acc_local78519[2])) / 2.2838315556293922983e-07);
    }
    {
      for (int slice_idx = 0; slice_idx < (((i78536 * 3) + 3) - (i78536 * 3)); slice_idx++) {
        sat_acc78457[((i78536 * 3) + slice_idx)] = (sat_acc78457[((i78536 * 3) + slice_idx)] + (2.8247609439046209905e-07 * grav_acc78520[(0 + slice_idx)]));
      }
    }
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        central_acc78458[(0 + slice_idx)] = (central_acc78458[(0 + slice_idx)] - (sat_gms78455[i78536] * grav_acc78520[(0 + slice_idx)]));
      }
    }
  }
  for (int i78671 = 0; i78671 < 4; i78671++) {
    {
      for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
        acc78521[(0 + slice_idx)] = (sat_acc78457[((i78671 * 3) + slice_idx)] - central_acc78458[(0 + slice_idx)]);
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i78671 * 6) + 3) - (i78671 * 6)); slice_idx++) {
        jupsatsystem_noderiv78449[((i78671 * 6) + slice_idx)] = state78451[(((i78671 * 6) + 3) + slice_idx)];
      }
    }
    {
      for (int slice_idx = 0; slice_idx < (((i78671 * 6) + 6) - ((i78671 * 6) + 3)); slice_idx++) {
        jupsatsystem_noderiv78449[(((i78671 * 6) + 3) + slice_idx)] = acc78521[(0 + slice_idx)];
      }
    }
  }
  return 0;
}

