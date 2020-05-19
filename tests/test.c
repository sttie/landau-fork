#include <math.h>

static inline long double get_dfdx_cell(int full_idx10691, int dx_mapped_size10692, int al_idx10693, int inv_mapping_period10694, int *restrict dx_idx_mappings10695, long double *restrict der_vec10696) {
  return ((al_idx10693 < inv_mapping_period10694) ? ((dx_idx_mappings10695[((inv_mapping_period10694 * full_idx10691) + al_idx10693)] >= 0) ? der_vec10696[((dx_mapped_size10692 * full_idx10691) + dx_idx_mappings10695[((inv_mapping_period10694 * full_idx10691) + al_idx10693)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx10691, int dx_mapped_size10692, int al_idx10693, int inv_mapping_period10694, int *restrict dx_idx_mappings10695, long double *restrict der_vec10696) {
  return ((al_idx10693 < inv_mapping_period10694) ? ((al_idx10693 == full_idx10691) ? 1.0 : ((dx_idx_mappings10695[((inv_mapping_period10694 * full_idx10691) + al_idx10693)] >= 0) ? der_vec10696[((dx_mapped_size10692 * full_idx10691) + dx_idx_mappings10695[((inv_mapping_period10694 * full_idx10691) + al_idx10693)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx10698, int inv_mapping_period10699, int *restrict dx_idx_mappings10700, long double *restrict der_vec10701) {
  return ((((al_idx10698 < inv_mapping_period10699) ? dx_idx_mappings10700[al_idx10698] : -1) >= 0) ? der_vec10701[((al_idx10698 < inv_mapping_period10699) ? dx_idx_mappings10700[al_idx10698] : -1)] : 0.0);
}

int y(long double *restrict y10702, long double m10703) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      y10702[(0 + slice_idx)] = m10703;
    }
  }
  return 0;
}

int f(long double *restrict f10708, long double *restrict m10709) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      f10708[(0 + slice_idx)] = m10709[(0 + slice_idx)];
    }
  }
  return 0;
}

int z(long double *restrict z10714, long double *restrict m10715, long double *restrict n10716) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      z10714[(0 + slice_idx)] = m10715[(0 + slice_idx)];
    }
  }
  return 0;
}

int g(long double *restrict g10721) {
  long double arr10723[3] = { 0.0 };
  {
    long double y10755[3] = { 0.0 };
    y(y10755, 1.0e0);
    long double f10762[3] = { 0.0 };
    f(f10762, arr10723);
    long double z10763[3] = { 0.0 };
    z(z10763, y10755, f10762);
    long double f10764[3] = { 0.0 };
    f(f10764, z10763);
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      arr10723[(0 + slice_idx)] = f10764[slice_idx];
    }
  }
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      g10721[(0 + slice_idx)] = arr10723[(0 + slice_idx)];
    }
  }
  return 0;
}

