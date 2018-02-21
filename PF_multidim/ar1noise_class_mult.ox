/*
ASIR inherited structure for AR(1) + noise model:
mkpitt 27/2/99
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include <oxprob.h>
#include "asir.ox"



// NEIL -- THIS IS NOT USED
// R from M;
// weight 1 * M
// output iy  1*R
deter_samp(const weight, const R)
{
  decl iy = zeros(1,R); 
  decl M = columns(weight);
  decl count = zeros(1,M);
  count = floor(R * weight);
  decl d = sumr(count);
  decl i ;
  decl next = 0.* ones(1,count[][0]);
  iy = next;
  for (i=1; i < M; i++) 
    {
      next = i.* ones(1,count[][i]);
      iy = iy ~ next;
    }
  decl Remainder = R -d;
  decl w_star = weight - count ./R ; 
  w_star =w_star./sumr(w_star);
  
  next = weighted_bootstrap(w_star, Remainder );
  //print(d~columns(iy)~R);
   //print(iy);
  //print(next); exit(1);
  //print(Remainder/R);
  return iy~next;
 
   
   
}


// Kalman filter: outputs amean, avar (filt mean and var for alpha_t | Y_t):
kf_filter(const amY,  const amean, const avar, const nobs, const phi, const q, 
	      const sig_eps, const a0, const p0, const mu)
{ 
/* Also adapted to return the log-likelihood */
  decl i=0; 
  decl a_star = a0, p_star = p0;
  decl a, p;
  decl l = zeros(1, nobs); 
  for (i=0; i < nobs; i++)
    {
       // alpha_i | Y_{i-1] ~ N( a*, p*)
       p = 1/(1/p_star + 1/sig_eps);
       a = p * ( amY[0][][i]/sig_eps + a_star/p_star);		// alpha_i | Y_{i-1] ~ N(
	   // alpha_i | Y_i ~ N( a, p)
       l[][i] = -0.5 * log( sig_eps + p_star ) - 0.5  *  log(2 * M_PI)
	      -0.5 *  (amY[0][][i]  - a_star)^2 / (sig_eps + p_star);
	   
       amean[0][][i] = a ;
       avar[0][][i] = p ;
       a_star = mu * (1 - phi) + phi * a;   // alpha_{i+1} | Y_i
       p_star = phi.^2 * p + q;	  	 
    }
return l;
}
// Kalman filter: outputs  (filt mean and var for alpha_t | Y_t):
kf2_filter(const mY,  const nobs, const phi, const q, 
	      const sig_eps, const a0, const p0, const mu)
{ 
/* Also adapted to return the log-likelihood */
  decl i=0; 
  decl a_star = a0, p_star = p0;
  decl a, p;
  decl l = zeros(1, nobs); 
  for (i=0; i < nobs; i++)
    {
       // alpha_i | Y_{i-1] ~ N( a*, p*)
       p = 1/(1/p_star + 1/sig_eps);
       a = p * ( mY[][i]/sig_eps + a_star/p_star);		// alpha_i | Y_{i-1] ~ N(
	   // alpha_i | Y_i ~ N( a, p)
//	   print(i~a_star~p_star);
       l[][i] = -0.5 * log( sig_eps + p_star ) //- 0.5  *  log(2 * M_PI)
	      -0.5 *  (mY[][i]  - a_star)^2 / (sig_eps + p_star);
	   
      
       a_star = mu * (1 - phi) + phi * a;   // alpha_{i+1} | Y_i
       p_star = phi.^2 * p + q;	  	 
    }
return l;
}

// Kalman filter: outputs  (filt mean and var for alpha_t | Y_t):
// model mY = mAlpha + N(0, mEPS2)
//       mAlpha = mT * mAlpha +  N(0, mHH); 
kfmult_filter(const mY,  const nobs, const mT, const mHH, 
	      const mEPS2, const mA0, const mP0, const mu)
{ 
/* Also adapted to return the log-likelihood */
  decl i=0; 
  decl mA_star = mA0, mP_star = mP0;
  decl mA, mP;
  decl dim = rows(mT);
  decl l = zeros(1, nobs); 
  for (i=0; i < nobs; i++)
    {
//	   print(i~mA_star~mP_star);
	   l[][i] = -0.5 * log( determinant(mP_star + mEPS2) ) 
	      -0.5 *  (mY[][i]  - mA_star)' * invertsym( mP_star + mEPS2) * (mY[][i]  - mA_star);
       // alpha_i | Y_{i-1] ~ N( a*, p*)
       mP = mP_star - mP_star * invertsym( mP_star + mEPS2) * mP_star ;   //1/(1/p_star + 1/sig_eps);
       mA = mA_star + mP_star * invertsym( mP_star + mEPS2) * (mY[][i] - mA_star);// a = p * ( mY[][i]/sig_eps + a_star/p_star);		// alpha_i | Y_{i-1] ~ N(
	   // alpha_i | Y_i ~ N( a, p)
//       l[][i] = -0.5 * log( sig_eps + p_star ) - 0.5  *  log(2 * M_PI)
//	      -0.5 *  (mY[][i]  - a_star)^2 / (sig_eps + p_star);
	   	
      
       mA_star = (unit(dim) - mT) * mu + mT * mA;   // alpha_{i+1} | Y_i
       mP_star = mT * mP * mT' + mHH;	  	 
    }
return l;
}



/* Self explan:	  */
Sim_AR1(const amY,  const amAlpha, const nobs, const mT, const sig_eta2, const sig_eps2,
        const mIP, const mu)
{
    decl i;
    decl dim = rows(amY[0]);
    amAlpha[0][][0] = 	mu + choleski(mIP) * rann(dim,1);//
    amY[0][][0] = amAlpha[0][][0] + sqrt(sig_eps2) * rann(dim,1);

    for (i=1; i<nobs; i++) {
		amAlpha[0][][i] = mu * (unit(dim) - mT)	  + mT * amAlpha[0][][i-1] +
		                  sqrt(sig_eta2) .* rann(dim,1);			
		//amY[0][][i] = exp(amAlpha[0][][i]/2) * rann(1,1);
		amY[0][][i] = amAlpha[0][][i]	+ sqrt(sig_eps2) * rann(dim,1);
	}
//	print(mT); print(mu~sig_eta2~ sig_eps2 ); exit(1);

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
        mExp[][t] = meanr(mAlpha_t) ; 
        quan[][t] = 0; //quantilec(exp(mAlpha_t./2)', <0.1, 0.3, 0.5, 0.7, 0.9>);
	    mVar[][t] =  varr(mAlpha_t);
        mlogf[][t] = logft;
	    return 0;
}




/* class inherited from Pfilter_Data */ 
//_______________________________________________________________START OF INHERITED SPEC______________________

class Pfilter_ar1_Data : Pfilter_Data
{
	 Pfilter_ar1_Data(const the_parameters, const s);
	 expand_pts(const mAlpha_t, const lag);
	 innov_func(const prev_state_t);
	 lik_func(const alpha_t,const meas_t);
	 
	 innov_lag1_func(const prev_state_t, const meas_t); // full adaption
	 lik_lag1_func(const mAlpha_t, const meas_tplus1);	// full adaption
	 innov_func_U(const mAprev_state, const mU, const mT)	;
	 initial_func(const M)	;
         predn_distance( const filter_store, const meas_t);
         sim_y(const mAlpha_prev, const amY_sim); 
	 get_params();
};	
// constructor. 
// In:   3 * 1 matrix of parameters.
Pfilter_ar1_Data::Pfilter_ar1_Data(const the_parameters, const the_s)
{
	 parameters = the_parameters;s = the_s;
}

// Simultes from density of initial state
// In:  integer of number of sampes
// returns:	1 * M matrix of simulations.
Pfilter_ar1_Data::initial_func(const M)
{
	decl mIP;
    if (parameters[1][] >= 1.0)  mIP = parameters[2][].^2; 
    else  mIP= parameters[2][]^2 /  (unit(s) - parameters[1][].^2) ;	
        
	decl vIa =  parameters[0][] ;
	return  vIa + choleski(mIP) * rann(s,M);
}
// Simulates 1 step-ahead
// In:   r * c matrix containing prev states
// Returns: r * c matrix containing c simulations of the states (each dim r)
// mU uniforms s* M.
Pfilter_ar1_Data::innov_func_U(const mAprev_state, const mU, const mT)
{
	decl M = columns(mAprev_state) ;
	return parameters[0][] + mT * (mAprev_state - parameters[0][])
	  + parameters[2][] .* mU;  
	
}
//Pfilter_ar1_Data::innov_func(const mAprev_state)
//{
//	decl M = columns(mAprev_state) ;
//	return parameters[0][]/s + parameters[1][] * (mAprev_state - parameters[0][]/s)
//	  + parameters[2][]/sqrt(s) .* rann(s,M); //quann(mU);  //(-log(ranu(1,M)) -1) ;	//	
//}


Pfilter_ar1_Data::lik_func(const mAlpha_t, const meas_t)
{
  decl M = columns(mAlpha_t),dim = rows(mAlpha_t);
  decl  answer = sumc(  ( meas_t - mAlpha_t ).^2 )  ./ parameters[3][].^2 ;
  return -0.5 * answer  -  0.5 * dim * log(  parameters[3][].^2 ) .* ones(1,M) ;	  //2 * M_PI *
}


Pfilter_ar1_Data::innov_lag1_func(const mAprev_state, const meas_t)
{
   decl M = columns(mAprev_state) ;
	// print(M);
   decl mMean =  parameters[0][] + parameters[1][] * (mAprev_state - parameters[0][]);
   decl sig_eta2 = parameters[2][]^2;
   decl sig_eps2 = parameters[3][].^2; //meas

   decl sigma2 = 1/(1/sig_eta2	+ 1/sig_eps2);
   decl mMean_post = sigma2 .* (meas_t/sig_eps2 + mMean ./ sig_eta2); 
   return  mMean_post + sqrt(sigma2) .* rann(1,M);	 

}
// for full adaption
Pfilter_ar1_Data::lik_lag1_func(const mAlpha_t, const meas_tplus1)
{
  decl M = columns(mAlpha_t);
  decl dim = rows(mAlpha_t);
  decl mMean = parameters[0][]/dim + parameters[1][] * (mAlpha_t - parameters[0][]/dim);
  decl ymean = sumc(mMean);
  
  decl yvar	= parameters[2][].^2 + parameters[3][].^2 ; // scalar..
  
  decl  answer = (meas_tplus1 - ymean).^2 ./ (yvar); // (1.837877066)
  return	  -0.5 * answer  -  0.5  *  log(2 * M_PI * (yvar) ) .* ones(1,M) ;
 //  decl answer = -0.5 .* mAlpha_t -0.5 .* meas_t^2 .* exp(-mAlpha_t);
  
  // return answer; 
}
/*
// Sim y(t) from prev alpha:
Pfilter_ar1_Data::sim_y(const mAlpha_prev, const amY_sim) 
{
    decl p = rows(amY_sim[0]);
    decl s = rows(mAlpha_prev);
    decl M = columns(mAlpha_prev);
    decl SIM = columns(amY_sim[0]);
    decl iy = ranu(1, SIM) * M;
    decl mAlphat = innov_func(mAlpha_prev[][iy]) ;
    //decl ft = exp(mAlphat[0][]/2) .* rann(1 ,SIM ); // 1*M
    //decl wt =  exp(mAlphat[1:(s-1)][]/2) .* rann(p ,SIM ); // p*M
    amY_sim[0][][] = mAlphat ; // CHANGE;
    return 0;
}
*/
/* Choose your expansion points.
// mAlpha_t: s * M  matrix of previous state points
// lag:  integer. How far ahead (0 = 1-step)  */
Pfilter_ar1_Data::expand_pts(const mAlpha_t, const lag)
{
	
        decl mAnswer;
// parameters[1] ^(lag+1)
	 mAnswer =  parameters[0][]+  parameters[1][] * (mAlpha_t - parameters[0][]);      
	return mAnswer;
}
Pfilter_ar1_Data::get_params()
{
	return parameters;
}
//_______________________________________________________________END OF INHERITED SPEC______________________


/*	
simulates Model 1:	 
time dim T, meas dim N, state dim k

mSige  	  in	N * N   matrix	 of   meas var
		  out 	unchanged
mZ		  in    N * k matrix
		  out	unchanged
mLambda	  in    k * k diagonal matrix
		  out   unchanged
mH		  in	k * k matrix of var Sign
		  out   unchanged
mIa		  in	k * 1 starting mean for states
		  out   unchanged

mIP		  in	k * k starting var for states
		  out   unchanged

amAlpha	  in	k * T matrix of zeros
		  out   k * T matrix of states
amSignal  in	N * T matrix of zeros
		  out   N * T matrix of signals
amY		  in    N * T matrix of zeros
		  out   N * T matrix of y's

*/


/*


SimMSV1( const mB, const mLambda, const mH,
		const mu, const mIP,
	    const amAlpha, const amY)
{
          decl N = 1; //rows(mB); //5
	  decl k = 1; //rows(mH); //6
	  decl T = columns(amY[0]);
	  
	  amY[0] = zeros(N, T);
	  amAlpha[0] = zeros(k, T);

	  decl mAlpha0 =  mu + choleski(mIP) * rann(k,1);
	  amAlpha[0][][0] = mu + mLambda * (mAlpha0-mu)  + mH * rann(k,1);;
          decl t = 0 ; 
          decl ft = exp(amAlpha[0][0][t]/2) * rann(1,1);
	  amY[0][][t] = mB * ft ;//+  diag(exp(amAlpha[0][1:(k-1)][t]/2))  * rann(N,1) ;	 
	  for (t=1; t < T; t++)
	  {
		  amAlpha[0][][t] = mu + mLambda * (amAlpha[0][][t-1 ] - mu) +
							mH * rann(k,1);
		   ft = exp(amAlpha[0][0][t]/2) * rann(1,1);
		  amY[0][][t] = mB * ft ;//+  diag(exp(amAlpha[0][1:(k-1)][t]/2))  * rann(N,1)  ;

	  }
	 
}

w0_sv(const mY, const mExpand_block, const t
		 , const filter_data, const amu_star      )
{
   decl expand = mExpand_block[0];
   decl M = columns(expand);
   decl sigma = filter_data->get_params()[2][];
   //print("sig", sigma);
   decl meas_t = mY[][t];
   decl exp_expand = (meas_t*meas_t) .* exp(-expand) ;
   decl mean_st = expand + 0.5 * sigma^2 * ( exp_expand   - ones(1, M));
   
   decl con1 = (0.5/sigma^2) .* (mean_st .* mean_st - expand .* expand) 
                - 0.5 * exp_expand  .* (ones(1, M) + expand);
   amu_star[0] = mean_st;
  
   decl lmax = 1.0 * max(con1);
  return exp(con1 - lmax);
  
}

w1_sv(const mY, const mAlpha, const mExpand_block, const t
		 , const filter_data)
{
  decl expand = mExpand_block;
  decl R = columns(mAlpha);
  decl meas_t = mY[][t]; //print(meas_t);
  
  decl weights = - 0.5 * (meas_t *meas_t ) .* ( exp(-mAlpha) - exp(-expand) .* 
      (ones(1, R)- (mAlpha - expand)) );
  //if (weights[][0] < -0.5) //printweights[][0]~t~meas_t~mAlpha[][0]~expand[][0]);
  decl lmax = 0.0; //1.0 * max(weights);
  //print(mAlpha'~expand'~weights');
  
  return exp(weights - lmax);
}
  
innov_adapt(const amSample_block, const iy, const filter_data, const mu_star)
{
   //decl R = columns(amSample_block[0][0]) ;
  decl R = columns(iy); //print(R);
   //print(filter_data->get_params());
   decl sigma = filter_data->get_params()[2][];// print(sigma);
   //decl mu =  filter_data->get_params()[0];
   decl answer = mu_star[][iy] +  sigma * rann(1, R);
   amSample_block[0][0]= answer;
   //print(columns(amSample_block[0][0])~R);
   //amSample_block[0][0][][] = mu_star[][iy] +  sigma * rann(1, R); 

}



asir_detersv(const mY,   const filter_data, const output_data, 
     const M, const SIR_OR_ASIR, const iFixedLag, const Y_SIM)
{
	decl T = columns(mY);	// time dim
	decl p = rows(mY);	 // obs dim
	decl s = filter_data->give_s(); // state dim	
	decl mAlpha = zeros(s, M);
	decl weights0 = zeros(1, M);
	decl weights1 = zeros(1,  M );
	decl i = 0;
	decl iy;
	decl t; 
	decl mExpand_pts;
	decl act_lik;
	decl expand_lik;
	decl filter_store = new array[iFixedLag+1];
	mAlpha = filter_data->initial_func(M);
	for (i=0; i<iFixedLag; i++)
	{
   		filter_store[i] = zeros(s, M);
	}
	filter_store[iFixedLag] = mAlpha;
	decl mSample_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mSample_block[i] = zeros(s,  M );
	}
	decl mExpand_block = new array[iFixedLag+1];
	for (i=0; i<=iFixedLag; i++)
	{
   		mExpand_block[i] = zeros(s, M);
	}

	decl lag = 0;
        decl ut = 0;
	decl loglt =0;
	decl loglt_sim = zeros(1, Y_SIM);
        decl mY_sim = zeros(p, Y_SIM); decl Uni; 
        decl mean_lt, varlt;  decl mu_star = zeros(s, M);
        decl ok = 0, it; it = 0;
        decl pi = ones(1,M);
	for (t =0; t < T; t++)
	{
          it = 0;
	 
	  if (t > 0) lag = min(t-1, iFixedLag);
	  expand_block((&mExpand_block), filter_store[iFixedLag - lag], lag, iFixedLag, filter_data);
	  weights0 = w0_sv(mY,mExpand_block,  t,  filter_data,  &mu_star);
          mExpand_block[0] = mu_star;
          weights0 = w0_sv(mY,mExpand_block,  t,  filter_data,  &mu_star);
	  weights0 = weights0 ./ sumr(weights0)	;
	  iy = deter_samp(weights0, M);
          
          //print(M~columns(iy));
	  // iy = weighted_bootstrap(weights0, M );   //p(weights0[][0:10], iy[][0:(20)], iy[][M-20:M-1]); 
          //print(weights0); 
	  // print(iy); 
          //exit(1);
	  innov_adapt(&mSample_block, iy, filter_data, mu_star);
          //loglt = predn_decomp(filter_store[iFixedLag],  mY[][t], filter_data, 1 ) ; 
          //print(t);
	  weights1 = w1_sv(mY,  mSample_block[0], mExpand_block[0][][iy],
				 t, filter_data);
          weights1 = weights1 ./ sumr(weights1) ; 
          iy = deter_samp(weights1, M);
	  mAlpha = mSample_block[iFixedLag]; 
	  mAlpha = mAlpha[][iy]; 
	  
      	  output_data->get_output(mAlpha, t, loglt);   
          mAlpha = mAlpha[][0:(M-1)];   
          // p(mAlpha[][0:10]);      
	  for (i=0; i<iFixedLag; i++)	filter_store[i] = filter_store[i+1];
	  filter_store[iFixedLag] = mAlpha;       // a_t | Y_t
         
       }

	delete filter_store, mExpand_block, mSample_block; 
	return 0;
 }

*/

