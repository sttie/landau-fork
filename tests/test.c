#include <math.h>

static inline long double get_dfdx_cell(int full_idx10728, int dx_mapped_size10729, int al_idx10730, int inv_mapping_period10731, int *restrict dx_idx_mappings10732, long double *restrict der_vec10733) {
  return ((al_idx10730 < inv_mapping_period10731) ? ((dx_idx_mappings10732[((inv_mapping_period10731 * full_idx10728) + al_idx10730)] >= 0) ? der_vec10733[((dx_mapped_size10729 * full_idx10728) + dx_idx_mappings10732[((inv_mapping_period10731 * full_idx10728) + al_idx10730)])] : 0.0) : 0.0);
}

static inline long double get_dfdx_cell_dx(int full_idx10728, int dx_mapped_size10729, int al_idx10730, int inv_mapping_period10731, int *restrict dx_idx_mappings10732, long double *restrict der_vec10733) {
  return ((al_idx10730 < inv_mapping_period10731) ? ((al_idx10730 == full_idx10728) ? 1.0 : ((dx_idx_mappings10732[((inv_mapping_period10731 * full_idx10728) + al_idx10730)] >= 0) ? der_vec10733[((dx_mapped_size10729 * full_idx10728) + dx_idx_mappings10732[((inv_mapping_period10731 * full_idx10728) + al_idx10730)])] : 0.0)) : 0.0);
}

static inline long double get_dfdx_var(int al_idx10735, int inv_mapping_period10736, int *restrict dx_idx_mappings10737, long double *restrict der_vec10738) {
  return ((((al_idx10735 < inv_mapping_period10736) ? dx_idx_mappings10737[al_idx10735] : -1) >= 0) ? der_vec10738[((al_idx10735 < inv_mapping_period10736) ? dx_idx_mappings10737[al_idx10735] : -1)] : 0.0);
}

int y(long double *restrict y10739, long double m10740) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      y10739[(0 + slice_idx)] = m10740;
    }
  }
  return 0;
}

int f(long double *restrict f10745, long double *restrict m10746) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      f10745[(0 + slice_idx)] = m10746[(0 + slice_idx)];
    }
  }
  return 0;
}

int z(long double *restrict z10751, long double *restrict m10752, long double *restrict n10753) {
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      z10751[(0 + slice_idx)] = m10752[(0 + slice_idx)];
    }
  }
  return 0;
}

int g(long double *restrict g10758) {
  long double arr10760[3] = { 0.0 };
  {
    long double y10792[3] = { 0.0 };
    y(y10792, 1.0e0);
    long double f10799[3] = { 0.0 };
    f(f10799, arr10760);
    long double z10800[3] = { 0.0 };
    z(z10800, y10792, f10799);
    long double f10801[3] = { 0.0 };
    f(f10801, z10800);
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      arr10760[(0 + slice_idx)] = f10801[slice_idx];
    }
  }
  {
    for (int slice_idx = 0; slice_idx < 3; slice_idx++) {
      g10758[(0 + slice_idx)] = arr10760[(0 + slice_idx)];
    }
  }
  return 0;
}

