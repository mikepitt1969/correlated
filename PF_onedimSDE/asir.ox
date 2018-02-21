	
/* ASIR.ox 31/07/97 Pitt & Shephard.
   General code to implement SIR/ASIR procedure
*/

#include <oxstd.h>
#include <oxdraw.h>

//#include "wboot.h"
#include "wboot.ox"
//decl g_jj;

 

// Class for ouputting results of filter:
class Pfilter_output
{
  	decl T	; // time dim
  	virtual get_output(const mAlpha_t, const t, const logft); // gives output 
};

// Class consisting of mainly virtual members allowing functions to be defined by
// Offspring. 
class Pfilter_Data
{
	 decl parameters;
	 decl s; // state dim
	 virtual innov_func(const mAlpha_t_prev); 
	 virtual lik_func(const mAlpha_t, const measurement_t);
	 virtual initial_func(const sampsz);
	 virtual expand_pts(const mAlpha_t, const lag);
     virtual sim_y(const mAlpha_prev, const amY_sim); 
	 give_s();
	 give_parameters();
	 reset_parameters(const the_parameters);
         
};
Pfilter_Data::give_parameters()  {return parameters;}
Pfilter_Data::give_s()  {return s;}
// Reset parameters: 
Pfilter_Data::reset_parameters(const the_parameters)
{
	parameters = the_parameters;
}

/*
// INPUT:
// mAlpha_prev
// mAlpha_prev, s * M matrix from previous time step.
// 
predn_decomp(const mAlpha_prev, const vY, const filter_data, const it)
{
    decl l = 0.0;
    decl mAlpha_t ;
    decl i = 0, lmean=0.0;
    if (it ==1)
    {
		mAlpha_t = filter_data->innov_func(mAlpha_prev);
		l = exp( filter_data->lik_func(mAlpha_t, vY));
        lmean = meanr(l);
    }
    else 
      { 
	for (i=0; i < it; i++)
	  {
	    mAlpha_t = filter_data->innov_func(mAlpha_prev);
	    l = exp( filter_data->lik_func(mAlpha_t, vY));
	    lmean += meanr(l);
	  }
	lmean /= it;
      }
    return log(lmean)';
}  
*/	 
// evaluation of likelihood:
lik_lag(const mY, const mAlpha_period, const time, const lag, const iFixed_lag, const filter_data)
{
	decl l = 0.0;
	decl i;
	for (i=(iFixed_lag - lag); i <= iFixed_lag; i++)
		l += filter_data->lik_func(mAlpha_period[i], mY[][time - iFixed_lag +i]);
	return exp(l);		 // - lmax
}
// choose expansion spray of points:
expand_block(const amExpand_block, const mAlpha_start, const lag,
			 const iFixed_lag, const filter_data)
{
	decl i =0;
	for (i=(iFixed_lag - lag); i <= iFixed_lag; i++)
	  amExpand_block[0][i] = filter_data->expand_pts(mAlpha_start, i - (iFixed_lag - lag) );
	return 0;

}
// simulate from innov equation ahead into mSample_block:
innov_lag(const amSample_block, const mAlpha_start, const lag,
		  const iFixedLag, const filter_data)
{
	decl i =0;
	decl mAlpha_prev = mAlpha_start	;
	for (i=(iFixedLag - lag); i <= iFixedLag; i++)
	{
		amSample_block[0][i] = filter_data->innov_func(mAlpha_prev);
		mAlpha_prev = amSample_block[0][i]	;
	}

	return 0;
}

/*
(Auxiliary)	Sampling Importance Resampling Filter
_________________________________________________

dimension of state -- s
dimension of meas  -- p
dimension of time  -- T


VARIABLE		INPUT								OUTPUT

mY			p * T matrix of observations		unchanged

filter_data		ptr to object of class inherited 	unchanged
				from Pfilter_Data containing 
				parameters and functions from 
				virtual parents
output_data		ptr to object of class inherited 	unchanged
				from Pfilter__output containing 
				parameters and functions from 
				virtual parents
M               int of number of simulations		unchanged
				propagated
R               int of number of simulations		unchanged
				resampled
SIR_OR_ASIR     int (1 for ASIR, 0 for SIR)         unchanged
iFixedLag	    int number of lags (0 for usual)    unchanged
LIK             int (1 for calculation of likelihhod, 0 for usual)


*/
asir(const mY,   const filter_data, const output_data, 
     const M, const R, const SIR_OR_ASIR, const iFixedLag, const LIK)
{
	decl T = columns(mY);	// time dim
	decl p = rows(mY);	 // obs dim
	decl s = filter_data->give_s(); // state dim	        
	decl mAlpha = zeros(s, M);
	decl weights0 = zeros(1, M);
	decl weights1 = zeros(1, R);
	decl i = 0;
	decl iy;
	decl t; 
	decl mExpand_pts;
	decl act_lik;
	decl expand_lik;
	decl filter_store = new array[iFixedLag+1];
	mAlpha = filter_data->initial_func(M);
	mAlpha=sortr(mAlpha);
	
	for (i=0; i<iFixedLag; i++)
	{
   		filter_store[i] = zeros(s, M);
	}
	filter_store[iFixedLag] = mAlpha;
	decl mSample_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mSample_block[i] = zeros(s, R);
	}
	decl mExpand_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mExpand_block[i] = zeros(s, M);
	}

	decl lag = 0;
        decl ut = 0;
	decl loglt =0;
    decl mean_lt, varlt;
	decl prob0, prob1; 
/* This is where all the action takes place */
	for (t =0;	t < T; t++)
	{
		if (t > 0) lag = min(t-1, iFixedLag);
		expand_block((&mExpand_block), filter_store[iFixedLag - lag], lag, iFixedLag, filter_data);	
               
		if (SIR_OR_ASIR == 0) prob0 = weights0 = 1/M .* ones(1, M);
		else {
			weights0 = lik_lag(mY, mExpand_block, t, lag, iFixedLag, filter_data);
			prob0 = weights0 ./ sumr(weights0);
		}                
		if  (SIR_OR_ASIR == 0) iy = ranu(1, R) .* M;
		else iy = weighted_bootstrap(prob0,R);              
       	innov_lag(&mSample_block, filter_store[iFixedLag - lag][][iy], lag, 
			  iFixedLag, filter_data);   			
		act_lik = lik_lag( mY, mSample_block, t, lag, iFixedLag, filter_data);
		if (SIR_OR_ASIR == 1) expand_lik = weights0[][iy]; 
		else  expand_lik = ones(1, R); 
		weights1 = act_lik ./ expand_lik ;
	/* THIS EVALUATES THE LIKELIHOOD */
		if (LIK == 1)
		{
			if (SIR_OR_ASIR == 0)  loglt = log(meanr(weights1)) ;
			else loglt = log(meanr(weights0)) +
			          log( meanr(weights1) ) + 0.5/R * varr(weights1)/meanr(weights1)^2; 
		}
		prob1 = weights1 ./ sumr(weights1) ;	 
		iy = weighted_bootstrap(prob1, M); 
		mAlpha = mSample_block[iFixedLag]; 
		mAlpha = mAlpha[][iy];
	   // Sort stage
		mAlpha=sortr(mAlpha);
              
	    output_data->get_output(mAlpha, t, loglt);       
		for (i=0; i<iFixedLag; i++)	filter_store[i] = filter_store[i+1];
		filter_store[iFixedLag] = mAlpha;       // a_t | Y_t
         
	}
	delete filter_store, mExpand_block, mSample_block; 
	return 0;
 }


/*
(Auxiliary)	Sampling Importance Resampling Filter
_________________________________________________

dimension of state -- s
dimension of meas  -- p
dimension of time  -- T


VARIABLE		INPUT								OUTPUT

mY			p * T matrix of observations		unchanged

filter_data		ptr to object of class inherited 	unchanged
				from Pfilter_Data containing 
				parameters and functions from 
				virtual parents
output_data		ptr to object of class inherited 	unchanged
				from Pfilter__output containing 
				parameters and functions from 
				virtual parents
M               int of number of simulations		unchanged
				propagated
R               int of number of simulations		unchanged
				resampled
SIR_OR_ASIR     int (1 for ASIR, 0 for SIR)         unchanged
iFixedLag		int number of lags (0 for usual)    unchanged
LIK             int (1 for calculation of likelihhod, 0 for usual)


*/
 
asir_strat(const mY,   const filter_data, const output_data, 
     const M, const R, const SIR_OR_ASIR, const iFixedLag, const LIK)
{
	decl T = columns(mY);	// time dim
	decl p = rows(mY);	 // obs dim
	decl s = filter_data->give_s(); // state dim	        
	decl mAlpha = zeros(s, M);
	decl weights0 = zeros(1, M);
	decl weights1 = zeros(1, R);
	decl i = 0;
	decl iy;
	decl t; 
	decl mExpand_pts;
	decl act_lik;
	decl expand_lik;
	decl filter_store = new array[iFixedLag+1];
	mAlpha = filter_data->initial_func(M);
	mAlpha=sortr(mAlpha);
	
	for (i=0; i<iFixedLag; i++)
	{
   		filter_store[i] = zeros(s, M);
	}
	filter_store[iFixedLag] = mAlpha;
	decl mSample_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mSample_block[i] = zeros(s, R);
	}
	decl mExpand_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mExpand_block[i] = zeros(s, M);
	}

	decl lag = 0;
        decl ut = 0;
	decl loglt =0;
    decl mean_lt, varlt;
	decl prob0, prob1; 	 decl  msort;
/* This is where all the action takes place */
	for (t =0;	t < T; t++)
	{
		if (t > 0) lag = min(t-1, iFixedLag);
		expand_block((&mExpand_block), filter_store[iFixedLag - lag], lag, iFixedLag, filter_data);	
               
		if (SIR_OR_ASIR == 0) prob0 = weights0 = 1/M .* ones(1, M);
		else {
			weights0 = lik_lag(mY, mExpand_block, t, lag, iFixedLag, filter_data);
			prob0 = weights0 ./ sumr(weights0);
		}                
		if  (SIR_OR_ASIR == 0) iy = weighted_strat_bootstrap(prob0,R);
		else iy = weighted_strat_bootstrap(prob0,R);              
       	innov_lag(&mSample_block, filter_store[iFixedLag - lag][][iy], lag, 
			  iFixedLag, filter_data);   			
		act_lik = lik_lag( mY, mSample_block, t, lag, iFixedLag, filter_data);
		if (SIR_OR_ASIR == 1) expand_lik = weights0[][iy]; 
		else  expand_lik = ones(1, R); 
		weights1 = act_lik ./ expand_lik ;
		
	/* THIS EVALUATES THE LIKELIHOOD */
		if (LIK == 1)
		{
			if (SIR_OR_ASIR == 0)  loglt = log(meanr(weights1)) ;
			else
				{
					loglt = log(meanr(weights0)) +
			                log( meanr(weights1) )+ 0.5/R * varr(weights1)/meanr(weights1)^2;; 
				}		
		}		
		mAlpha = mSample_block[iFixedLag]; // 1 * R
		prob1 = weights1 ./ sumr(weights1) ; // 1* R
		
		msort = sortbyr(mAlpha | prob1, 0);
		mAlpha = 	msort[0][];
		prob1 = msort[1][]; 
		
		iy = weighted_strat_bootstrap(prob1, M); 	  
		mAlpha = mAlpha[][iy];
		
		mAlpha=sortr(mAlpha);
   	    output_data->get_output(mAlpha, t, loglt);       
		for (i=0; i<iFixedLag; i++)	filter_store[i] = filter_store[i+1];
		filter_store[iFixedLag] = mAlpha;       // a_t | Y_t
               
	}
	delete filter_store, mExpand_block, mSample_block; 
	return 0;
 }


// 1* R, 1 * (M), 1*R, 1 * (M+1).
adjust_explik(const mSample, const mghat, const mindex, const mcells)
{
  decl num =  ( mcells[][mindex+1] - mSample ) .* mghat[][mindex]
           +  ( mSample - mcells[][mindex] ) .* mghat[][mindex+1] ;
  
  decl den =  ( mcells[][mindex+1] - mcells[][mindex] );
			  
  return (num ./ den); 

}
// 	Perfectly simulate from dens prop f(y|a) * f(a|a_{i-1})
perfect_sim(const yt, const mAlpha, const filter_data)
{
    decl M = columns(mAlpha);
	decl mSample = zeros(1, M), a, b;
	decl mU = zeros(1,M), logl, i, mean, Phi_a, Phi_b, trunc_u, trunc_n;
	mSample = mAlpha;
	decl parameters =  filter_data->give_parameters();
	decl mSE = parameters[3][], mSEta = parameters[2][],
	     mu = parameters[0][], mT =parameters[1][];
	decl constant = 0.5  *  log(2 * M_PI * (mSE.^2) ) .* ones(1,M);
	for (i=0; i < 4; i++)
	{
	     logl = filter_data->lik_func(mSample, yt)
		      + constant ;
		 // [u | x]
	   	 mU = exp(logl) .* ranu(1,M) ;
		 
		 // [x | u]
		 a = yt - mSE * sqrt(-2 * log(mU));
		 b = yt + mSE * sqrt(-2 * log(mU));	// 1* M
		 mean = mu + mT * (mAlpha - mu)	; // 1* M
		 // Check
		 
		 //print(log(mU[][0:5]) | filter_data->lik_func(a[][0:5], yt) + constant[][0:5]);
		 //print(log(mU[][0:5]) | filter_data->lik_func(b[][0:5], yt) + constant[][0:5]);
		 Phi_a = probn((a - mean) ./ mSEta );
		 Phi_b = probn((b - mean) ./ mSEta );
		// if (i > 3) print(Phi_a[][M-10:M-1] | Phi_b[][M-10:M-1] | mean[][M-10:M-1] | mSE);
		 trunc_u = ranu(1,M) .* (Phi_b - Phi_a) +  Phi_a; // 1* M
		 trunc_n = quann(trunc_u);
		 mSample = mean + mSEta * trunc_n ;
		 
		 
	}
	return mSample;
}
// ASIR:
Smasir1_strat(const mY,   const filter_data, const output_data, 
     const M, const R,  const iFixedLag, const LIK)
{
	decl T = columns(mY), p = rows(mY),  s = filter_data->give_s(); // state dim	        
	decl mAlpha = zeros(s, M),  weights0 = zeros(1, M),  weights1 = zeros(1, R);
	decl i = 0,  iy,  t,  mExpand_pts, act_lik, expand_lik, lag = 0, loglt =0 ;
	decl prob0, prob1,   msort,  mSample, mSim;
	
	decl filter_store = new array[iFixedLag+1];	 
	for (i=0; i<iFixedLag; i++)  filter_store[i] = zeros(s, M);	 
	decl mSample_block = new array[iFixedLag+1], mExpand_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mSample_block[i] = zeros(s, R);
		mExpand_block[i] = zeros(s, M); 
	}	
	mAlpha = filter_data->initial_func(M);
	mAlpha=sortr(mAlpha);
	filter_store[iFixedLag] = mAlpha;
	decl h = 0.01 * sqrt(1000)/sqrt(R), smooth_prob1; 
/* This is where all the action takes place */
	for (t =0;	t < T; t++)
	{
		if (t > 0) lag = min(t-1, iFixedLag);
		expand_block((&mExpand_block), filter_store[iFixedLag - lag], lag, iFixedLag, filter_data);	       
		weights0 = lik_lag(mY, mExpand_block, t, lag, iFixedLag, filter_data);
		prob0 = weights0 ./ sumr(weights0);
	//	iy = weighted_strat_bootstrap(prob0,R);
	    mSample = weighted_smstrat_bootstrap(prob0,mAlpha,R);
		mAlpha = mSample[1][];
		iy = mSample[0][];
		
		innov_lag(&mSample_block, mAlpha, lag, iFixedLag, filter_data);
	//	innov_lag(&mSample_block, filter_store[iFixedLag - lag][][iy],  lag, iFixedLag, filter_data);
		mAlpha = mSample_block[iFixedLag]; 

		
		act_lik = lik_lag( mY, mSample_block, t, lag, iFixedLag, filter_data);
		weights1 = act_lik ./ weights0[][iy] ;
		msort = sortbyr(mAlpha | weights1, 0);
		mAlpha = msort[0][];	//t	not t+1
		weights1 = msort[1][];
	//	print(weights1[][0:10]);
		weights1 =  sumr(weights1) * pre_smooth(mAlpha, h, weights1);
		//print(weights1[][0:10]);
		//print(smooth_prob1[][0:10]);
		if (LIK == 1) loglt = 	//	log(meanr(weights0)) + log( meanr(weights1) );
		                  //     + 0.5/R * varr(weights1)/meanr(weights1)^2;
		//log( meanr( lik_lag(mY, mSample_block, t, lag, iFixedLag, filter_data)));



	    /* 
		prob1 = weights1 ./ sumr(weights1) ; // 1* R	
		msort = sortbyr(mAlpha | prob1, 0);
		mAlpha = msort[0][];	//t	not t+1
		prob1 = msort[1][];
	
		smooth_prob1 =  pre_smooth(mAlpha, h, prob1);//	prob1; */
		smooth_prob1 = weights1 / sumr(weights1);
	 	mAlpha = weighted_smstrat_bootstrap(smooth_prob1,mAlpha,M)[1][];	   	
   	    output_data->get_output(mAlpha, t, loglt);	 		
		for (i=0; i<iFixedLag; i++)	filter_store[i] = filter_store[i+1];
		filter_store[iFixedLag] = mAlpha;       // a_t | Y_t		
	//	print(weights0[][0:20]);
	//	print(weights1[][0:20]);
	//	print(prob1[][0:20]*R);		
	}
	delete filter_store, mExpand_block, mSample_block; 
	return 0;
 }

// SIR 
Smasir0_strat(const mY,   const filter_data, const output_data, 
     const M, const R,const iFixedLag, const LIK)
{
	decl T = columns(mY),  p = rows(mY),  s = filter_data->give_s(); // state dim	        
	decl mAlpha = zeros(s, M),  weights0 = zeros(1, M),  weights1 = zeros(1, R);
	decl i = 0,  iy,  t,  act_lik,  expand_lik;
	decl lag = 0,  ut = 0, loglt =0, prob0, prob1, msort,  mSample,  mSim;
	decl filter_store = new array[iFixedLag+1];
	decl mSample_block = new array[iFixedLag+1];
	for (i=0; i<iFixedLag; i++)  filter_store[i] = zeros(s, M);	  
	for (i=0; i<=iFixedLag; i++)  mSample_block[i] = zeros(s, R);	
	
	mAlpha = filter_data->initial_func(M);
	mAlpha=sortr(mAlpha);
	filter_store[iFixedLag] = mAlpha;
/* This is where all the action takes place */
	for (t =0;	t < T; t++)
	{
	// SET SEEDS
	//    ranseed(g_jj +t);  setmyseed(19 + g_jj + t);
	// SEEDS SET
		if (t > 0) lag = min(t-1, iFixedLag);		      
		prob0 = weights0 = 1/M .* ones(1, M);
		iy = weighted_strat_bootstrap(prob0,R);
		mAlpha = filter_store[iFixedLag - lag][][iy];
		// Alt
		//mAlpha = weighted_smstrat_bootstrap(prob0,mAlpha,R)[1][]; 
		
       	innov_lag(&mSample_block, mAlpha, lag, 
			  	iFixedLag, filter_data);
		act_lik = lik_lag( mY, mSample_block, t, lag, iFixedLag, filter_data);
	
		weights1 = act_lik  ; 
		if (LIK == 1)	loglt = log(meanr(weights1)) + 0.5/R * varr(weights1)/meanr(weights1)^2;;	 
		mAlpha = mSample_block[iFixedLag]; // t+1, 1 * R
		prob1 = weights1 ./ sumr(weights1) ; // 1* R		
		msort = sortbyr(mAlpha | prob1, 0);
		mAlpha = 	msort[0][];
		prob1 = msort[1][];
	 	mAlpha = weighted_smstrat_bootstrap(prob1,mAlpha,M)[1][];
   	    output_data->get_output(mAlpha, t, loglt);	 
		for (i=0; i<iFixedLag; i++)	filter_store[i] = filter_store[i+1];
		filter_store[iFixedLag] = mAlpha;       // a_t | Y_t
	}
	delete filter_store, mSample_block; 
	return 0;
 }

/* NEW - NOV 07	- simplified functions for methods 1 to 4
 Mtd 1 - 1 step GSS
 Mtd 2 - 1 step perfect adaption
 Mtd 3 - 1 step global partial adaption
 Mtd 4 - local partial adaption ASIR  
*/

/* Nov 07 - GSS type mtd - Does not involve R */ 
PF_loop2(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, weights; 
   // time 0:
   mAlpha = filt_obj->initial_func(M); mAlpha = sortr(mAlpha);	 
   for (t= 0; t < T; t++)
   {
    	mAlpha = filt_obj->innov_func(mAlpha); mAlpha = sortr(mAlpha);	// 1* M	
		weights = exp(filt_obj->lik_func(mAlpha, mY[][t]));
		loglt = log(meanr(weights))+ 0.5/M * varr(weights)/meanr(weights)^2;
		weights = weights ./ sumr(weights) ; // 1* R
		mAlpha = weighted_smstrat_bootstrap(weights,mAlpha,M)[1][];
		out_obj->get_output(mAlpha, t, loglt);	
   }
   return 0; 
}


/* Nov 07 - full adapt mtd - Does not involve R */ 
PF_loop_fulladapt(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{

//   lik_lag1_func(const mAlpha_t, const meas_tplus1)
// innov_lag1_func(const mAprev_state, const meas_t)
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, weights; 
   // time 0:
   mAlpha = filt_obj->initial_func(M); mAlpha = sortr(mAlpha);	 
   for (t= 0; t < T; t++)
   {
		weights = exp(filt_obj->lik_lag1_func(mAlpha, mY[][t]));  
		loglt = log(meanr(weights)); //+ 0.5/M * varr(weights)/meanr(weights)^2;
		weights = weights ./ sumr(weights) ; // 1* R
//		mAlpha = weighted_smstrat_bootstrap(weights,mAlpha,M)[1][];

		// DELETE THIS NEXT BIT - WE WANT SMOOTH
		if (SM_OR_NOT  == 0)   iy = weighted_bootstrap(weights, M);
		else   iy = weighted_strat_bootstrap(weights,M); 

		mAlpha = mAlpha[][iy];
		//________________________________
		mAlpha = filt_obj->innov_lag1_func(mAlpha, mY[][t]); mAlpha = sortr(mAlpha);	// 1* M	
		out_obj->get_output(mAlpha, t, loglt);	
   }
   return 0; 
}

/* Nov 07 - global approximate adapt mtd - Does not involve R */ 
PF_loop_globaladapt(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{

//   lik_lag1_func(const mAlpha_t, const meas_tplus1)
// innov_lag1_func(const mAprev_state, const meas_t)
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, loglt2 = 0.0, weights, weights0, expand;
   decl mAlpha_inter,mAlpha_new;
   // time 0:												 
   mAlpha = filt_obj->initial_func(M); mAlpha = sortr(mAlpha);
   mAlpha_inter	= mAlpha;
   decl probs = 1/M .* ones(1,M); 
   for (t= 0; t < T; t++)
   {
//		mAlpha_new = filt_obj->innov_func(mAlpha_inter);  mAlpha_new = sortr(mAlpha_new);	
//		weights = exp(filt_obj->lik_func(mAlpha_new, mY[][t]));
//		loglt2 = log(meanr(weights) );
//		weights = weights ./ sumr(weights) ; // 1* R
//		mAlpha_new = weighted_smstrat_bootstrap(weights,mAlpha_new,M)[1][];




		
		expand = meanr(mAlpha); // expansion pt, change
		weights = exp(filt_obj->appr_lik_lag1_func(mAlpha, mY[][t], expand));	
		weights = weights .* (probs);  // not prperly weighted here
		loglt = log(sumr(weights)); 
//		print("mtd now ", loglt);
//		weights = exp(filt_obj->appr_lik_lag1_func(mAlpha_inter, mY[][t], expand));
//		loglt = log(meanr(weights)); 
//		print("mtd 2", loglt);
//		
//		mAlpha_new = filt_obj->innov_func(mAlpha_inter);
//		weights = exp(filt_obj->appr_lik_func(mAlpha_new, mY[][t], expand));
//		loglt = log(meanr(weights)); 
//		print("mtd 0", loglt);	
		
		weights = weights ./ sumr(weights) ; // 1* M
		mAlpha = weighted_smstrat_bootstrap(weights,mAlpha,M)[1][];
		mAlpha = filt_obj->approx_innov_lag1_func(mAlpha, mY[][t], expand);
		mAlpha = sortr(mAlpha);	// 1* M	
		probs = exp(filt_obj->lik_func( mAlpha, mY[][t])
		        - filt_obj->appr_lik_func(mAlpha, mY[][t], expand));
		loglt += log(meanr(probs)) + 0.5/M * varr(probs)/meanr(probs)^2;; 
		probs =  probs ./ sumr(probs) ; // 1* M
		mAlpha_inter = weighted_smstrat_bootstrap(probs,mAlpha,M)[1][];	 //filtered	 t | t
//		print("malpha mtd 1 ", sumr(mAlpha .* probs));
		out_obj->get_output(mAlpha_inter, t, loglt);	//  _inter
//		print(loglt~loglt2);
//		exit(1);
   }
//   exit(1);
   return 0; 
}
	 



/* Nov 07 - global approximate adapt mtd - Does not involve R */ 
PF_loop_globaladapt2(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{

//   lik_lag1_func(const mAlpha_t, const meas_tplus1)
// innov_lag1_func(const mAprev_state, const meas_t)
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, loglt2 = 0.0, weights, weights0, expand;
   decl mAlpha_inter;
   // time 0:
   mAlpha = filt_obj->initial_func(M); mAlpha = sortr(mAlpha);
   mAlpha_inter	= mAlpha;
   decl probs = 1/M .* ones(1,M); 
   for (t= 0; t < T; t++)
   {
//   		weights = exp(filt_obj->lik_func(mAlpha_inter, mY[][t]));
//		expand = meanr(mAlpha_inter); // expansion pt
//        weights0 = exp(filt_obj->appr_lik_lag1_func(mAlpha_inter, mY[][t], expand));
//		loglt = log(meanr(weights0)); //+ 0.5/M * varr(weights0)/meanr(weights0)^2;
		loglt2 = log(meanr(	exp(filt_obj->lik_func(mAlpha_inter, mY[][t])) ) );
		print(filt_obj->lik_func(mAlpha_inter, mY[][t]) );
		expand = meanr(mAlpha); // expansion pt
		weights = exp(filt_obj->appr_lik_lag1_func(mAlpha, mY[][t], expand));
		
		weights = weights .* (probs);  // not prperly weighted here
		loglt = log(sumr(weights)); 
//		loglt = log(meanr(weights))+ 0.5/M * varr(weights)/meanr(weights)^2;
		weights = weights ./ sumr(weights) ; // 1* M
		mAlpha = weighted_smstrat_bootstrap(weights,mAlpha,M)[1][];
		mAlpha = filt_obj->approx_innov_lag1_func(mAlpha, mY[][t], expand);
		mAlpha = sortr(mAlpha);	// 1* M	
		probs = exp(filt_obj->lik_func( mAlpha, mY[][t])
		        - filt_obj->appr_lik_func(mAlpha, mY[][t], expand));
		loglt = log(meanr(probs)); //+ 0.5/M * varr(probs)/meanr(probs)^2;
		probs =  probs ./ sumr(probs) ; // 1* M
		print(filt_obj->lik_func(mAlpha, mY[][t]));
		mAlpha_inter = weighted_smstrat_bootstrap(probs,mAlpha,M)[1][];
		// f(alpha{t+1} | Y{t+1} )
//		print(meanr(mAlpha_inter)~meanr(mAlpha));
	   // or pass weights as well! 
		out_obj->get_output(mAlpha_inter, t, loglt);	//  _inter
		print(loglt~loglt2);
   }
   exit(1);
   return 0; 
}
	 