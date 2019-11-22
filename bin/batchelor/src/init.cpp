#include "batchelor.h"

#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    // MNN calculations.
    REGISTER(find_mutual_nns, 2),
    REGISTER(cosine_norm, 2),
    REGISTER(smooth_gaussian_kernel, 4),
    REGISTER(adjust_shift_variance, 6),

    // Unlogging calculations.
    REGISTER(unlog_exprs_mean, 4),
    REGISTER(unlog_exprs_scaled, 5),
    {NULL, NULL, 0}
};

void attribute_visible R_init_batchelor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

