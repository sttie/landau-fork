#include <math.h>

static inline long double get_dfdx_cell(int full_idx8957, int dx_mapped_size8958, int al_idx8959, int inv_mapping_period8960, int *restrict dx_idx_mappings8961, long double *restrict der_vec8962) {
  return ((al_idx8959 < inv_mapping_period8960) ? ((dx_idx_mappings8961[((inv_mapping_period8960 * full_idx8957) + al_idx8959)] >= 0) ? der_vec8962[((dx_mapped_size8958 * full_idx8957) + dx_idx_mappings8961[((inv_mapping_period8960 * full_idx8957) + al_idx8959)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx8957, int dx_mapped_size8958, int al_idx8959, int inv_mapping_period8960, int *restrict dx_idx_mappings8961, long double *restrict der_vec8962) {
  return ((al_idx8959 < inv_mapping_period8960) ? ((al_idx8959 == full_idx8957) ? 1.0 : ((dx_idx_mappings8961[((inv_mapping_period8960 * full_idx8957) + al_idx8959)] >= 0) ? der_vec8962[((dx_mapped_size8958 * full_idx8957) + dx_idx_mappings8961[((inv_mapping_period8960 * full_idx8957) + al_idx8959)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx8964, int inv_mapping_period8965, int *restrict dx_idx_mappings8966, long double *restrict der_vec8967) {
  return ((((al_idx8964 < inv_mapping_period8965) ? dx_idx_mappings8966[al_idx8964] : -1) >= 0) ? der_vec8967[((al_idx8964 < inv_mapping_period8965) ? dx_idx_mappings8966[al_idx8964] : -1)] : 0.0);
}

int f(long double *restrict f8968, long double *restrict m8969) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      f8968[(0 + slice_idx)] = 1.0e0;
    }
  }
  return 0;
}

int g(long double g8974) {
  long double f8976;
  f8976 = 1.0e0;
  g8974 = f8976;
  return 0;
}

