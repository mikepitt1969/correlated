#include <oxstd.h>


extern "asirsmooth,C_multinomial_ran"
Ranmultinomial(const n, const weight, x);

extern "asirsmooth,C_strat_multinomial_ran"
Ranstratmultinomial(const n, const weight, x, u);

extern "asirsmooth,C_smooth1"
Ransmooth1(const iy, const mX, const mu, const weights, mY);

extern "asirsmooth,C_PreSmooth"
PreSmooth(const mX, const weight, h,  smooth_weight); //const h,

extern "asirsmooth,C_smooth2"
Ransmooth2(const iy, const mX, const mu, const weights, mY);

extern "asirsmooth,C_settheseed"
setmyseed(const my_seed);

extern "asirsmooth,C_strat_multinomial_ran2"
Ranstratmultinomial2(const n, const weight, const u1, x);

 