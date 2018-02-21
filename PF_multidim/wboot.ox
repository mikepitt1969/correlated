#include <oxstd.h>
#include "wboot.h"

// INPUT:
// weight -- a  1 * M vector (giving normalised multinomial probabilities i.e. +ve , sum to 1).
// i_number  -- an integer of required number of samples from the multinomial.
//
pre_smooth(const mx, const h, const weight)	// , const h
{
   decl smooth_weight;
   PreSmooth(mx, weight,  h, &smooth_weight);// h,
   return smooth_weight ./ sumr(smooth_weight);
}


// RETURNS:
// A 1 * (i_number) vector of the sampled indices. 
weighted_bootstrap(const weight, const i_number)
{
	decl iy ;
	Ranmultinomial(i_number,weight,&iy);
	return iy;
}

// RETURNS:
// Hilbert Sorting
// A 1 * (i_number) vector of the sampled indices. 
hilbert_sort(const mX)
{
	decl iy ;
	Hilbert_Srt(mX, &iy);
	return iy;
}

// INPUT:
// weight -- a  1 * M vector (giving normalised multinomial probabilities i.e. +ve,
//				sum to 1).
// i_number  -- an integer of required number of samples from the multinomial.
//
// RETURNS:
// A 1 * (i_number) vector of the sampled indices. 
weighted_strat_bootstrap(const weight, const i_number)
{
	decl iy, mu ;
	Ranstratmultinomial(i_number,weight,&iy, &mu);
	return iy;
}
// new funtion - February 2015
// we pass u1 - a scalar Uniform variate and the only source of randomness.
weighted_strat_bootstrap2(const weight, const i_number, const u1)
{
	decl iy, mu ;
	Ranstratmultinomial2(i_number,weight,u1,&iy);
	return iy;
}
kernel(const vx,  const avwage, 
			  const h, const prob)
{
	decl answer, bit1, bit2, i;
	decl n = rows(avwage); 

	answer = prob .* densn((vx - avwage)/h) * 1/n;	  
	answer = sumc(answer);

	return answer;

}  
check_alpha(const mPr,const mX)
{
   decl R = columns(mX);
   decl d = 0.001;
   decl a = 15/d;
   decl mProb = mPr;	 
   decl mDiff = mX[][1:(R-1)] - mX[][0:(R-2)]; decl weight ;
   decl ind = vecindex(mDiff .<= d), i, index, diff, add	;
  // mProb = kernel(
   if (rows(ind) > 0)
   {  //print("index");
	  for (i=0; i < rows(ind); i++)
	  {
	      index = ind[i][];
		  
	      diff = mX[][index] - mX[][index - 1];
	      weight =   exp( a * (diff - d/2) );
		  weight = 1 / (1 + weight);
		  add = 0.5 * ( mPr[][index+1] - mProb[][index] ) * weight;
		  mProb[][index] = mProb[][index] + add;
		  mProb[][index +1 ] = mProb[][index+1] - add;
		 // print(index~weight~mPr[][index]~mProb[][index]);
		  
	  }

	}
   /*
   	 weight =   exp( a .* (mDiff[][ind] - d/2) );  
     weight = 1 ./ (1 + weight);  
    
    mProb[][ind] = mPr[][ind] + 0.5 .* ( mPr[][ind+1] - mPr[][ind] ) .* weight;
    mProb[][ind +1 ] = mPr[][ind+1] - 0.5 .* ( mPr[][ind+1] - mPr[][ind] ) .* weight;
   }  */
   if (sumr(mProb) < 0.9999999 && rows(ind) > 0) {
   //  print(sumr(mProb));
//	 print( ind' | mDiff[][ind] | mX[][ind] | mX[][ind+1] | weight);
	 
   }
   return mProb ;
   
}
// THIS DOES SMOOTH BOOTSRAPPING OF R dim vector
// INPUT:
// weight -- a  1 * R vector (giving normalised multinomial probabilities i.e. +ve,
//				sum to 1).
// i_number  -- an integer of required number of samples from the multinomial.
//
// RETURNS:
// A 1 * (i_number) vector of the sampled indices. 
weighted_smstrat_bootstrap(const weight, const mX, const i_number)
{
	decl iy, mu;
	decl mY ;
	decl R = columns(mX); 
	decl cat_weights = zeros(1, R-1);
	cat_weights[][0] = 0.5 * (2 * weight[][0] + weight[][1]);
	cat_weights[][R-2] = 0.5 * ( weight[][R-2] + 2 * weight[][R-1]);
	cat_weights[][1:(R-3)] = 0.5 * (weight[][1:(R-3)] + weight[][2:(R-2)]);
	cat_weights	=  cat_weights ./ sumr(cat_weights);
	Ranstratmultinomial(i_number,cat_weights,&iy, &mu);
	Ransmooth1(iy, mX, mu, weight,  &mY);
	return (iy | mY);
}


// THIS DOES SMOOTH BOOTSRAPPING OF R dim vector
// INPUT:
// expan -- a  1 * R vector (giving normalised multinomial probabilities i.e. +ve,
//				sum to 1).
// i_number  -- an integer of required number of samples from the multinomial.
//
// RETURNS:
// A 1 * (i_number) vector of the sampled indices. 
weighted_smquadstrat_bootstrap(const weight, const mX, const i_number)
{
	decl iy, mu;
	decl mY ;
	decl R = columns(mX); 
	decl cat_weights = zeros(1, R-1);
	cat_weights[][0] = 0.5 * (2 * weight[][0] + weight[][1]);
	cat_weights[][R-2] = 0.5 * ( weight[][R-2] + 2 * weight[][R-1]);
	cat_weights[][1:(R-3)] = 0.5 * (weight[][1:(R-3)] + weight[][2:(R-2)]);
	cat_weights = cat_weights ./ sumr(cat_weights);
	Ranstratmultinomial(i_number,cat_weights,&iy, &mu);
	Ransmooth2(iy, mX, mu, weight,  &mY);
	return (iy | mY);
}

