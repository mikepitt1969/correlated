/*
ASIR inherited structure for AR(1) + noise model:
mkpitt 27/2/99
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include <oxprob.h>
#include "asir.ox"
//	 C:\Users\Michael\my_files\PAPERS_PRIORITY2015\slow_Z_copula_paper\New_code_Nov15\SSF_Hilbert
// copied from ar1noise_class_pierre



/* Self explan:	vol  */
Sim_VAR1(const amY,  const amAlpha, const nobs, const mT, const sig_eta2, const rho,
        const mIP, const mu)
{
    decl i;
    decl dim = rows(amY[0]), innov;
    amAlpha[0][][0] = 	mu + choleski(mIP) * rann(dim,1);//
    amY[0][][0] =  exp(amAlpha[0][][0]/2)  * rann(dim,1);

    for (i=1; i<nobs; i++) {
	    innov =  rann(dim,1);
		amAlpha[0][][i] = mu * (unit(dim) - mT)	  + mT * amAlpha[0][][i-1] +
		                  sqrt(sig_eta2) .* innov;			
	
		amY[0][][i] = exp(amAlpha[0][][i]/2) * 
		      ( rho * innov + sqrt(1-rho^2) * rann(dim,1) );//	 * rann(dim,1) ;
	}
    return 0;
}

class Pfilter_ar1_output : Pfilter_output
{
	decl T	;
	decl mExp ;
	decl mVar ;decl mlogf,  mlogfnorm, m_u, mdiff, quan;
	decl s;
	Pfilter_ar1_output(const T, const s);
	get_output(const mAlpha_t, const t, const logft); 
		 //  const mlogfnorm, const m_u, const mdiff);
	give_mExp();
	give_mVar();give_mlogf()  ;give_mlogfn();give_mu(); give_mdiff();
        give_quan();
};
Pfilter_ar1_output::give_quan()	{return quan  ;}
Pfilter_ar1_output::give_mExp()	{return mExp  ;}
Pfilter_ar1_output::give_mVar()	{return mVar  ;}
Pfilter_ar1_output::give_mlogf()        {return mlogf;}
Pfilter_ar1_output::give_mlogfn() { return mlogfnorm ;}
Pfilter_ar1_output::give_mu() { return m_u;}
Pfilter_ar1_output::give_mdiff() { return mdiff;}
Pfilter_ar1_output::Pfilter_ar1_output(const the_T, const the_s)
{
	 T = the_T ;
	 s = the_s ;
	 mExp = zeros(s, T);
	 mVar = zeros(s, T);
     mlogf = zeros(1, T);
     mlogfnorm = zeros(1, T);
     m_u = zeros(1, T);
     mdiff = zeros(1, T);
     quan = zeros(5, T);
         
}
// CURRENT vYt, NEXT vYtn
// CAN BE ALTERED -- DEPENDS ON WHAT IS REQ'D:
Pfilter_ar1_output::get_output(const mAlpha_t, const t, const logft) 
			//	 const mlogftnorm, const ut, const diff )
{
        mExp[][t] = meanr(exp(mAlpha_t/2)) ; 
        quan[][t] = 0; //quantilec(exp(mAlpha_t./2)', <0.1, 0.3, 0.5, 0.7, 0.9>);
	    mVar[][t] =  varr(mAlpha_t);
        mlogf[][t] = logft;
	    return 0;
}




/* class inherited from Pfilter_Data */ 
//_______________________________________________________________START OF INHERITED SPEC______________________

class Pfilter_ar1_Data : Pfilter_Data
{
	 
	 decl Eul_delta	;
	 Pfilter_ar1_Data(const the_parameters, const s, const the_Eul_data);
	 expand_pts(const mAlpha_t, const lag);
	 innov_func(const prev_state_t);
	 lik_func(const alpha_t,const meas_t, const mSigma2_t, const mGamma_t);
	 
	 innov_lag1_func(const prev_state_t, const meas_t); // full adaption
	 lik_lag1_func(const mAlpha_t, const meas_tplus1);	// full adaption
	 innov_func_U(const mAprev_state, const mU)	;
	 give_Eul_delta();
	 initial_func(const M)	;
         predn_distance( const filter_store, const meas_t);
         sim_y(const mAlpha_prev, const amY_sim); 
	 get_params();
};	
// constructor. 
// In:   3 * 1 matrix of parameters.
Pfilter_ar1_Data::Pfilter_ar1_Data(const the_parameters, const the_s, const the_Eul_data)
{
	 parameters = the_parameters;s = the_s;
	 Eul_delta = the_Eul_data; 
}
 Pfilter_ar1_Data::give_Eul_delta() {return Eul_delta	;	}
// Simultes from density of initial state
// In:  integer of number of sampes
// returns:	1 * M matrix of simulations.

Pfilter_ar1_Data::initial_func(const mU)
{

   decl phi = parameters[1][];	decl sig = parameters[2][]	;
   decl mu = parameters[0][];
   decl nu = 2 * mu * (1-phi) /	sig^2;
   decl beta = 2 * (1-phi)/	sig^2;
   
   decl mUnif = probn(mU); //uniform
   decl mX = quangamma(mUnif, nu, beta); //print(mX);
   return log(mX);

}



// Now the CIR/SQRT model:
// from OU
mu_func(const mX, const parameters)
{
	decl phi = parameters[1][];	decl sig = parameters[2][]	;
    decl kappa = 1- phi ;
	decl mu = parameters[0][];
	decl answer = kappa .* ( mu-exp(mX)) .* exp(-mX) -  sig^2/2 .* exp(-mX);

	return  answer;  
}

sigma2_func(const mX, const parameters)
{								  
	decl sig = parameters[2][] ;
	decl answer = sig^2 .* exp(-mX);
//	answer[][vecindex(mX[0][] .< -10)] = 10e-6;
    return answer ;
}

 
// Simulates 1 step-ahead
// In:   1 * M matrix containing prev states
//		mU uniforms M_Eul * M.
// Returns: 1 * M matrix containing c simulations of the states (each dim r)
// 
Pfilter_ar1_Data::innov_func_U(const mAprev_state, const mU)
{
	decl M = columns(mAprev_state) ;
	decl mu = parameters[0][], phi = parameters[1][] , sig = parameters[2][] ; 
 	decl time_steps = rows(mU);
	decl mX = mAprev_state, v;
	decl SUM = exp(mX), SUM2 = 0, mXo;
	if (any(mX .< -10))
	{
		v = vecindex( mX[0][] .< -10);	
		mX[][v] = -10;
	    print("innov mX", mX);//print(	v );	print(mX); exit(1);
	}
	for (decl i=0; i < time_steps; i++)
	{
	    mXo = mX;
		mX = mX + mu_func(mX, parameters) .* Eul_delta
				+ sqrt(Eul_delta) .* sqrt(sigma2_func(mX, parameters)) .* mU[i][];
	    if (any(mX .< -10))	{ v = vecindex( mX[0][] .< -10);	mX[][v] = -10;	}
		if (isnan(mX)) {print("innov step", mXo, parameters, mu_func(mXo, parameters), sigma2_func(mXo, parameters)); exit(1); }
		
		SUM2 += exp(mX/2) .* mU[i][];
		if (i < (time_steps -1)) 	SUM += exp(mX);	  
	}	
	decl mSigma2 = Eul_delta * SUM;
	decl mGamma = sqrt(Eul_delta) * SUM2;
	
	return mX | mSigma2 | mGamma;
}



//innov_func(const mAprev_state)
Pfilter_ar1_Data::lik_func(const mAlpha_t, const meas_t, const mSigma2_t, const mGamma_t)
{
  decl rho = parameters[3][]; //leverage 
  decl M = columns(mAlpha_t);
//  decl mVar = (1-rho^2) .* exp(mAlpha_t)  , mInVar = 1/(1-rho^2) .* exp(-mAlpha_t);
  decl mMean = 0;  //correct this.
  // OR continuous time version
  decl answer = -0.5 * log(2 * M_PI) -0.5 * log(1-rho^2)-0.5 .* log(mSigma2_t)
        -0.5 .* (meas_t   - rho * mGamma_t).^2
			./ ( mSigma2_t * (1 - rho^2) );
  if (isnan(answer))
  {
  	print("like", answer,rho,mSigma2_t,mGamma_t);
	exit(1);
 }
//  decl  answer =  meas_t.^2 .* mInVar ;
//  answer = -0.5 * answer  -  0.5 * log( 2 * M_PI .* mVar ) .* ones(1,M);
  return answer ;
}

// parameters:
// mu, theta (0.98), sig_eta, rho
log_prior(const parameters)
{
 // 	parameters[0][0] = mu,  = ,		  //mu
  //	parameters[2][0] = sig_eta, parameters[3][0] = rho;
  decl mu =	parameters[0][0] , theta  =	parameters[1][0],
  sig_eta= parameters[2][0]  , rho= parameters[3][0]  ;
  return -0.5 * (mu-1)^2/0.2^2 -0.5 * (theta - 0.97)^2/0.005^2
    -0.5 * (log(sig_eta^2) - log(0.03))^2/0.1 -0.5 * (rho +0.65)^2/0.05^2;

}





Pfilter_ar1_Data::get_params()
{
	return parameters;
}

 Hsortr(const mX);
/* june 2016  - GSS type mtd - Does not involve R
 This is either stratified or not
 Involves U
 Now for the Euler scheme 
 mU	an array T+1 * M_Eul * (M+1)   // last element for the resample part
 */ 
PF_loop_U_sort(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT, const mU)
{
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl dim = s, u1;
   decl mAlpha = zeros(s,M), t = 0, iy=range(0, M-1,1), loglt = 0.0, weights;
   decl parameters= filt_obj->give_parameters();
   mAlpha =	filt_obj->initial_func(mU[0][0][1:M]);
   mAlpha = sortr(mAlpha);
   
   decl mFilter_t, mSigma2_t=0,  mGamma_t=0; 
   for (t= 0; t < T; t++)
   {

//       iy = Hsortr(mAlpha);//range(0, M-1,1);
//		mAlpha = mAlpha[][iy];		
//		print(t~ mAlpha[0][0:10])	; 
//		print(T~t~mAlpha)	 ;	
    	mFilter_t = filt_obj->innov_func_U(mAlpha,mU[t+1][][1:M]);		 	// 1 * M: Innovation state equation
		mAlpha	=  mFilter_t[0][] ;	mSigma2_t = mFilter_t[1][] ;	mGamma_t = 	mFilter_t[2][] ;

   	
//		iy = Hsortr(mAlpha);//range(0, M-1,1); //	//   	iy = range(0, M-1,1);
//			print(iy);
			
//		print(iy); exit(1);
//		if (t >= 3)  exit(1);

		iy = sortcindex(mAlpha);
		mAlpha = mAlpha[][iy]; mSigma2_t = mSigma2_t[][iy] ; mGamma_t = mGamma_t[][iy]	 ;
		weights = exp( filt_obj->lik_func(mAlpha, mY[][t], mSigma2_t,  mGamma_t) );
		loglt = log(meanr(weights));	   //		print(" var: ", varr(log(weights[0][])) ) ;
		weights = weights ./ sumr(weights) ; // 1* R
		u1 = probn(mU[t+1][0][0])	;  // part for resapmpling:	 //		print(weights, M, u1);
		iy = weighted_strat_bootstrap2(weights, M, u1);
		mAlpha = mAlpha[][iy];	 //		print(t~ iy)  ;
		out_obj->get_output(mAlpha, t, loglt);	
   }
 
   return 0; 
}

/* june 2016  - GSS type mtd - Does not involve R
 This is either stratified or not
 Involves U
 Now for the Euler scheme 
 mU	an array T+1 * M_Eul * (M+1)   // last element for the resample part
 */ 
PF_loop_U_nosort(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{
   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl dim = s, u1;
   decl mAlpha = zeros(s,M), t = 0, iy=range(0, M-1,1), loglt = 0.0, weights;
   decl parameters= filt_obj->give_parameters();
   decl Eul_delta= filt_obj->give_Eul_delta();
//   decl mu = parameters[0][], theta = parameters[1][], sig_eta2 = parameters[2][].^2,rho = parameters[3][];
   mAlpha =	filt_obj->initial_func(rann(1,M));
   //mu + sqrt(sig_eta2/(1-theta^2))  .* quann(mU[0][0][1:M]);	 // 1*M
//   print(T~mAlpha)	 ;
//   mAlpha = Hsortr(mAlpha);
   
   decl mFilter_t, mSigma2_t=0,  mGamma_t=0, mU; 
   for (t= 0; t < T; t++)
   {
		mU = rann(1/Eul_delta, M);
    	mFilter_t = filt_obj->innov_func_U(mAlpha,mU);		 	// 1* M
		mAlpha	=  mFilter_t[0][] ;	mSigma2_t = mFilter_t[1][] ;	mGamma_t = 	mFilter_t[2][] ;

		weights = exp( filt_obj->lik_func(mAlpha, mY[][t], mSigma2_t,  mGamma_t) );
		loglt = log(meanr(weights));	   //		print(" var: ", varr(log(weights[0][])) ) ;
		weights = weights ./ sumr(weights) ; // 1* R
		u1 = ranu(1,1)	;  // part for resapmpling:	 //		print(weights, M, u1);
		iy = weighted_strat_bootstrap2(weights, M, u1);
		mAlpha = mAlpha[][iy];	 //		print(t~ iy)  ;
		out_obj->get_output(mAlpha, t, loglt);	
   }
   return 0; 
}
