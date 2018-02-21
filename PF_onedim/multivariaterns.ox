
/* General algorithm for partially observed SDEs */
#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include <oxprob.h>
#include "multivariaterns.h"

ranmultT(const vmu, const mSig, const dof)
{
   decl m = rows(vmu);
   decl propose = vmu + choleski(mSig)
				* rann(m,1) ./ sqrt( rangamma(1,1,dof/2, dof/2) );
   return propose;

}
	 

/* standard (non-copular) multivariate t-density  */


logdensmultT(const x, const vmu, const mSig, const dof)
{
    decl p = rows(x), ok = 0;
	decl kernel = (x - vmu)' * invertsym(mSig) * (x - vmu);
	//if (invertsym(mSig) == 0) print(mSig);
	decl answer = -(dof + p)/2 * log(1 + kernel/dof)
			- 0.5 * logdet(mSig, &ok) ;
	return answer;
}