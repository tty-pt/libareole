/* Unity amalgamation entry point for libareole.
 * Preserves single-TU compilation compatibility across native, WASM, and unit tests
 * while keeping each internal subsystem decoupled in its own focused C unit. */

#include "areole_internal.h"
#include "ar_math.c"
#include "ar_occ.c"
#include "ar_pool.c"
#include "ar_populate.c"
#include "ar_solver.c"
#include "ar_probe.c"
