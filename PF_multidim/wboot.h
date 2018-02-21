#include <oxstd.h>


extern "HilbertCode,C_multinomial_ran"
Ranmultinomial(const n, const weight, x);

extern "HilbertCode,C_strat_multinomial_ran"
Ranstratmultinomial(const n, const weight, x, u);

extern "HilbertCode,C_smooth1"
Ransmooth1(const iy, const mX, const mu, const weights, mY);

extern "HilbertCode,C_PreSmooth"
PreSmooth(const mX, const weight, h,  smooth_weight); //const h,

extern "HilbertCode,C_smooth2"
Ransmooth2(const iy, const mX, const mu, const weights, mY);

extern "HilbertCode,C_settheseed"
setmyseed(const my_seed);

extern "HilbertCode,C_strat_multinomial_ran2"
Ranstratmultinomial2(const n, const weight, const u1, x);

extern "HilbertCode,C_HilbertSrt"
Hilbert_Srt(const x, index);





 