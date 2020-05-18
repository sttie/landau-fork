#include <math.h>

static inline long double get_dfdx_cell(int full_idx9259, int dx_mapped_size9260, int al_idx9261, int inv_mapping_period9262, int *restrict dx_idx_mappings9263, long double *restrict der_vec9264) {
  return ((al_idx9261 < inv_mapping_period9262) ? ((dx_idx_mappings9263[((inv_mapping_period9262 * full_idx9259) + al_idx9261)] >= 0) ? der_vec9264[((dx_mapped_size9260 * full_idx9259) + dx_idx_mappings9263[((inv_mapping_period9262 * full_idx9259) + al_idx9261)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx9259, int dx_mapped_size9260, int al_idx9261, int inv_mapping_period9262, int *restrict dx_idx_mappings9263, long double *restrict der_vec9264) {
  return ((al_idx9261 < inv_mapping_period9262) ? ((al_idx9261 == full_idx9259) ? 1.0 : ((dx_idx_mappings9263[((inv_mapping_period9262 * full_idx9259) + al_idx9261)] >= 0) ? der_vec9264[((dx_mapped_size9260 * full_idx9259) + dx_idx_mappings9263[((inv_mapping_period9262 * full_idx9259) + al_idx9261)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx9266, int inv_mapping_period9267, int *restrict dx_idx_mappings9268, long double *restrict der_vec9269) {
  return ((((al_idx9266 < inv_mapping_period9267) ? dx_idx_mappings9268[al_idx9266] : -1) >= 0) ? der_vec9269[((al_idx9266 < inv_mapping_period9267) ? dx_idx_mappings9268[al_idx9266] : -1)] : 0.0);
}

int f(long double *restrict f9270, long double *restrict m9271) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      f9270[(0 + slice_idx)] = 1.0e0;
    }
  }
  return 0;
}

int g(long double* g9276) {
  long double f9278;
  f9278 = 1.0e0;
  *g9276 = f9278;
  return 0;
}

