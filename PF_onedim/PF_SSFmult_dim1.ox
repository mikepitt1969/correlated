// 16/02/99 mkpitt. Full test of AR(1) +noise model via Monte Carlo. 
// For four scenarios. 
#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include "ar1noise_class_mult.ox"
#include "multivariaterns.ox"
//#include "kf_functions.ox"
//#include "multivariaterns.ox"
decl g_mY, g_mSig_eps, g_mL0, g_sig_eta;

PF_loop_U_nosort(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT);
set_mP0(const mT, const  sig_eta );
set_mP0(const mT, const  sig_eta );
Hsortr(const mX);
PF_loop_U(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT, const mU);
/* Do the AR(1) + noise examples */ 
get_returns(const data);
set_mT(const theta, const dim);
set_mP0(const mT, const  sig_eta ) ;
lik_kf(const vPar, const adFunc, const avScore, const amHessian);
exact_MCMC(const aExact_output, const mu_mean, const mu_var, const rho);
particle_MCMC(const aPart_output, const mu_mean, const mu_var,
		const rho, const N, const phi, const state_dim, const SIM_ini, const SIM_fix, const SIM_update,
		const filename);
/*
   Run exact MH (for SIM0 draws) and CPMH (for SIM draws) for dim = 1, theta = 0.4 for model of paper.
   T=400, dim = 1 and single parameter
   The MH and CPMH output is saved in two *.mat files for later analysis
   Values for beta = 0.47 and psi = 0.5 (see Table and Figure from paper)
*/
main()
{  
  decl theta = 0.4, dim = 1, T =1600, SIM = 20000, SIM0 = 10000;   //SIM0 updates for standard MH
  decl N , delta, phi,alpha = dim/(dim+1), SIM_ini = 100, SIM_fix = 2000, SIM_update =200; // initial (acc = 1), par fix, updates
  N = int(  0.47 * T^alpha ); delta =  0.5 * N/T; //0.12	
  phi = exp(-delta);  print("N~phi ", N~phi);
 
  decl mu = 0, sig_eta = 1, mSig_eps = sqrt(1);	 //sqrt(1-theta^2)
  decl mT = set_mT(theta, dim),  mP0=set_mP0( mT,  sig_eta ); 	
  decl mY = zeros(dim, T),  mAlpha_true = zeros(dim, T), mAlpha = mAlpha_true;	   
  ranseed(30); Sim_AR1( &mY,  &mAlpha_true, T, mT, sig_eta.^2, mSig_eps.^2, mP0, mu );
  
  decl parameters = ones(4,1);
  parameters[0][0] = mu, parameters[1][0] = theta, parameters[2][0] = sig_eta, parameters[3][0] = mSig_eps;
  g_mY = mY, g_mSig_eps = mSig_eps, g_sig_eta=sig_eta;	// MaxControl(160, 0);
  decl filt_obj = new Pfilter_ar1_Data(parameters, dim), out_obj = new Pfilter_ar1_output(T, dim) ;
 	 
  decl exact_output = zeros(1, SIM0), EX_IF = zeros(1, dim), EX_ACC= zeros(1, 1), rho = 0.7;
  decl particle_output= zeros(1+2, SIM); //contains theta and two more elements
  decl PF_IF = EX_IF, PF_ACC = EX_ACC ;	decl vPar = 0.413, mHess = 0.4 /T; 
  EX_ACC =  exact_MCMC(&exact_output, vPar, mHess, rho);
  savemat("mPar_MH1.mat",exact_output'~EX_ACC .* ones(SIM0, 1) ) ;
  PF_ACC = particle_MCMC(&particle_output, vPar, mHess, rho, N, phi, dim,
  		SIM_ini, SIM_fix, SIM_update, "mPar_CPM1.mat");
  Draw(0,  exact_output,"",  1, 1, 1);	Draw(1,  particle_output[0][],"",  1, 1, 1);
  /* output has been saved to the two files listed above - can then be analysed in another ox program */
  ShowDrawWindow();
  
}






// in mX is d * M
Hsortr(const mX)
{
//  decl m = meanr(mX);
//  decl V = sqrt(varr(mX));
//  decl mXs = (mX - m)./V ;
//  decl mU = exp( mXs)./(1+exp( mXs));
//  decl iy =  hilbert_sort(mU')	 ;
  return sortr(mX); 
}

set_mT(const theta, const dim)
{
  decl mT = zeros(dim, dim);
  for (decl i =0; i < dim; i++) { 
  	for (decl j=0; j < i; j++) mT[i][j] = mT[j][i] = theta^(i-j+1); //
	mT[i][i] = theta ;
  }
  return mT;
}
		

 

/* june 2016  - GSS type mtd - Does not involve R
 This is either stratified or not
 Involves U 
 mU	an array T+1 * dim * (M+1)
 */ 
PF_loop_U(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT, const mU)
{

   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl dim = s, u1;
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, weights;
   decl parameters= filt_obj->give_parameters();
   decl mu = parameters[0][], theta = parameters[1][], sig_eta2 = parameters[2][].^2,sig_eps2 = parameters[3][].^2;
   // time 0: mAlpha = filt_obj->initial_func(M);
   
   decl mT = set_mT(theta, dim);   decl mP0=set_mP0( mT,  parameters[2][] );  decl mL0= choleski(mP0);
   mAlpha =	mu + mL0  * (mU[0][][1:M]);	//  quann

//   print(mT,	 theta ); exit(1);
   mAlpha = Hsortr(mAlpha);	//Hsortr 
   for (t= 0; t < T; t++)
   {

    	mAlpha = filt_obj->innov_func_U(mAlpha,mU[t+1][][1:M], mT );		 	// 1* M
		mAlpha = Hsortr(mAlpha);  //Hsortr
		weights = exp(filt_obj->lik_func(mAlpha, mY[][t]));
		loglt = log(meanr(weights));	   //		print(" var: ", varr(log(weights[0][])) ) ;
		weights = weights ./ sumr(weights) ; // 1* R
		u1 = probn(mU[t+1][0][0])	;  // part for resapmpling:	 //		print(weights, M, u1);
		iy = weighted_strat_bootstrap2(weights, M, u1);
		mAlpha = mAlpha[][iy];	 //
//		if (t <10) print(t~ mAlpha[][0:9])  ;
//		else exit(1);
		out_obj->get_output(mAlpha, t, loglt);	
   }
   return 0; 
}

 
// simple state space model
/*
  y_t = I'alpha_t + sig_eps * eps_t
  alpha_t+1 = mu/d + phi * (alpha_t - mu/d) + sig/sqrt(d) * u_t
  */
 

 set_mP0(const mT, const  sig_eta ) 
{
	decl dim = rows(mT);
	decl mP0 = unit(dim) ; 
    for (decl i =0; i < 10; i++) { 
 	    mP0 = mT * mP0 * mT' + sig_eta.^2 .* unit(dim) ;
	}
//  exit(1);
  return mP0; 

}
// still to do ...
lik_kf(const vPar, const adFunc, const avScore, const amHessian)
{
   decl dim = rows(g_mY), T = columns(g_mY);;
   decl theta = vPar[0][0];//
   decl sig_eps = g_mSig_eps;
   decl mT =  set_mT(theta, dim); 
   decl mP0	= set_mP0(mT, g_sig_eta );   
   decl l =	kfmult_filter(g_mY,  T, mT, g_sig_eta^2 * unit(dim), 
	      g_mSig_eps^2 * unit(dim), zeros(dim,1), mP0, zeros(dim,1))   ;
//   kf2_filter(g_mY,  T,  mT, mH^2 , mSig_eps.^2, mu, mIP, mu);  
   adFunc[0] = sumr(l);
   return 1 ;   

}




/* june 2016  - GSS type mtd - Does not involve R
 This is either stratified or not
 Involves U 
 mU	an array T+1 * dim * (M+1)
 */ 
PF_loop_U_nosort(const filt_obj,  const out_obj, const mY,
		const M,  const SM_OR_NOT)
{

   decl T = columns(mY),  p = rows(mY),  s = filt_obj->give_s(); // state dim
   decl dim = s, u1;
   decl mAlpha = zeros(s,M), t = 0, iy, loglt = 0.0, weights;
   decl parameters= filt_obj->give_parameters();
   decl mu = parameters[0][], theta = parameters[1][], sig_eta2 = parameters[2][].^2,sig_eps2 = parameters[3][].^2;
   // time 0: mAlpha = filt_obj->initial_func(M);
   mAlpha =	mu + g_mL0 * rann(dim, M);	//  quann
   decl mT = set_mT(theta, dim); 

//   print(mT,	 theta ); exit(1);
//   mAlpha = Hsortr(mAlpha);	 
   for (t= 0; t < T; t++)
   {
//        mAlpha = Hsortr(mAlpha);
    	mAlpha = filt_obj->innov_func_U(mAlpha,rann(dim,M), mT );		 	//mU[t+1][][1:M] 1* M
//		mAlpha =  parameters[2][] .* (mU[t+1][][1:M]);	//quann
//		mAlpha = Hsortr(mAlpha);
		weights = exp(filt_obj->lik_func(mAlpha, mY[][t]));
		loglt = log(meanr(weights));	   //		print(" var: ", varr(log(weights[0][])) ) ;
		weights = weights ./ sumr(weights) ; // 1* R
		u1 = ranu(1,1)	;  // part for resapmpling:	 //		print(weights, M, u1);
		iy = weighted_strat_bootstrap2(weights, M, u1);
		mAlpha = mAlpha[][iy];	 //		print(t~ iy)  ;
		out_obj->get_output(mAlpha, t, loglt);	
   }
   return 0; 
}




exact_MCMC(const aExact_output, const mu_mean, const mu_var, const rho)
{	
	decl SIM = columns(aExact_output[0]), dim = rows(aExact_output[0]) ;
	decl ratio1, ratio2, curr_loglik = -10e6 ,new_loglik ;
	decl dF =0;
	decl vPar = 0.4 .*ones(dim, SIM), accept = 0.0, dof = 10;
	decl cond_mean1, cond_mean2, cond_var;
	
//	decl var = mu_var ;//* (dof - 2)/dof ;
	decl mL = choleski(mu_var);
	for (decl i=1 ; i < SIM; i++)
	{
		cond_mean1 = rho * vPar[][i-1] + (1-rho) * mu_mean	;
		vPar[][i] = cond_mean1 + sqrt(1-rho^2) * mL
				* rann(dim,1) ./ sqrt( rangamma(1,1,dof/2, dof/2) ); //vmu[0][]
		cond_mean2 = rho * vPar[][i] + (1-rho) * mu_mean	;
		cond_var = (1-rho^2) * mu_var ;
		
		ratio1 = logdensmultT(vPar[][i], cond_mean1, cond_var, dof)//	-0.5 * (vPar[][i] +0.4)^2/0.01; 
		       - logdensmultT(vPar[][i-1],cond_mean2 , cond_var , dof);// +0.5 * (vPar[][i] +0.4)^2/0.01 ;
		
		lik_kf(vPar[][i], &dF, 0, 0);
		new_loglik = dF;

		ratio2 = new_loglik - curr_loglik	;
		if (ranu(1,1) > exp(ratio2 - ratio1)) {
			vPar[][i] =	vPar[][i-1] ;//		
			accept += 1; //rejected
		}
		else  curr_loglik = new_loglik;
		if (fmod(i, 500) ==0 ) print("i exact: ", i, "\n");
	}
	aExact_output[0] = 	vPar ;
//	print("rej rate EXACT: ", 1- accept/SIM);
	return (1-accept/SIM);
}



//new_loglik = est_lik(vPar[0][i],  N);
//aPart_output[0] = 	vPar ;

particle_MCMC(const aPart_output, const mu_mean, const mu_var,
		const rho, const N, const phi, const state_dim, const SIM_ini, const SIM_fix, const SIM_update,
		const filename)
{	
	decl SIM = columns(aPart_output[0]), par_dim = rows(aPart_output[0]) -2;
	decl T = columns(g_mY), mL = choleski(mu_var);
	 
	decl ratio1, ratio2, curr_loglik = -10e6 ,new_loglik , dF =0;
	decl vPar = 0.4 * ones(par_dim, SIM), accept = 0.0, dof = 6;
	decl Z_prop = zeros(1, SIM), Z_acc = zeros(1, SIM), tru_loglik, est_loglik = zeros(1, SIM);
	decl vacc = zeros(1,SIM), vW = zeros(1,SIM), vWc = vW, W=0, W_curr=0;
	decl cond_mean1, cond_mean2, cond_var, parameters = ones(4,1);   //mu, mT, mH, mSig_eps
	decl mU = new array[T+1]; for (decl t=0; t < (T+1); t++) mU[t] =  rann(state_dim, N+1);
	decl mU_curr = new array[T+1];	for (decl t=0; t < (T+1); t++) mU_curr[t] = mU[t] ;

	parameters[0][0] = 0, parameters[1][0] = 0.4, parameters[2][0] = g_sig_eta,
  	parameters[3][0] = g_mSig_eps;
	
    decl filt_obj = new Pfilter_ar1_Data(parameters, state_dim);	
    decl out_obj = new Pfilter_ar1_output(T, state_dim) ;
 
	for (decl i=1 ; i < SIM; i++)
	{
		cond_mean1 =  rho * vPar[][i-1] + (1-rho) * mu_mean	;
		vPar[][i] =cond_mean1 + sqrt(1-rho^2) * mL
				* rann(par_dim,1) ./ sqrt( rangamma(1,1,dof/2, dof/2) ); //vmu[0][]
		
		cond_mean2 = rho * vPar[][i] + (1-rho) * mu_mean	;
		cond_var = (1-rho^2) * mu_var ;
		
		ratio1 = logdensmultT(vPar[][i], cond_mean1, cond_var, dof) //	+0.5 * (vPar[][i] +0.4)^2/0.01; 
		       - logdensmultT(vPar[][i-1],cond_mean2 , cond_var , dof);// +0.5 * (vPar[][i] +0.4)^2/0.01 ;
		if (i <= SIM_fix ) {vPar[][i] = 0.4; ratio1 = 0;  }
	    parameters[1][0] = 	vPar[][i] 	;
		
		filt_obj->reset_parameters(parameters);
		for (decl t=0; t < T; t++) mU[t] =  phi .* mU_curr[t]	 + sqrt(1-phi.^2) .* rann(state_dim,N+1) ;
		if (i <= SIM_ini ) {for (decl t=0; t < T; t++) mU[t] = rann(state_dim,N+1) ;} 
		PF_loop_U(filt_obj, out_obj, g_mY,	N,  1, mU)	;
		new_loglik = sumr( out_obj->give_mlogf() );

		ratio2 = new_loglik - curr_loglik	;
		
		lik_kf(vPar[][i], &dF, 0, 0);	tru_loglik = dF;
//		print(i~new_loglik~ tru_loglik~parameters')	;
		Z_prop[][i] = new_loglik - tru_loglik;
//		print(i~vPar[][i]~vPar[][i-1]~new_loglik~curr_loglik~Z_prop[][i]  )	;
//		print(ratio1~ratio2);
        W = Z_prop[][i] - Z_acc[][i-1] ;
		
		if (i <= SIM_ini ) {ratio2=ratio1;} 
		if (ranu(1,1) > exp(ratio2 - ratio1)) {
			vPar[][i] =	vPar[][i-1] ;//
			Z_acc[][i] = Z_acc[][i-1];
			accept += 1; //rejected
		}
		else //acccpt
		{
			curr_loglik = new_loglik;
			Z_acc[][i] = Z_prop[][i];
			W_curr = W;
//			mU_curr = mU;
			for (decl t=0; t < T; t++) mU_curr[t] =  mU[t] ;
		}
		if ( i > 2) vacc[][i] = 1-accept/(i+1) ;
		est_loglik[][i] = curr_loglik;
		vW[][i] = W;  vWc[][i] = W_curr;
		if (fmod(i, SIM_update) ==0 ) {print("acc rate PF: ",i~ 1-accept/(i+1)~T~N);
		      Draw(0, vPar[][0:i], 0, 1);
			  Draw(1, Z_prop[][0:i] | Z_acc[][0:i], 0, 1);
			  Draw(2, est_loglik[][1:i], 0, 1);
			  Draw(3, vacc[][1:i], 0, 1);
			  Draw(4, vW[][1:i] | vWc[][1:i], 0, 1);
			  ShowDrawWindow();
			  savemat(filename, vPar[][0:i]'~Z_prop[][0:i]'~Z_acc[][0:i]'~est_loglik[][0:i]'
			  	~vacc[][0:i]'~vW[][0:i]'~vWc[][0:i]') ;
			  }
//		print(i~ 1-accept/(i+1));
	}
	aPart_output[0] = 	vPar | Z_prop | Z_acc ;

//	print("rej rate PF: ", accept/SIM);

	return (1-accept/SIM);
}









 /*	 if (i > 0 && fmod(i, 20000) ==0 ) print("i, N, T, phi ",i~N~T~phi);  //  kf_lik = kf_filter(&mY, &mKF_mean, &mKF_var, T,  mT, sig_eta^2 , mSig_eps.^2, mu, mIP, mu);
//  print("kf_lik", kf_lik);

 if (i > 0 && fmod(i, 20000) ==0 ) print("i, N, T, phi ",i~N~T~phi);
	//	  kf_lik = kf_filter(&mY, &mKF_mean, &mKF_var, T,  mT, sig_eta^2 ,  mSig_eps.^2, mu, mIP, mu);
	//	  parameters[0][0] = mu;   filt_obj->reset_parameters(parameters);	 print(filt_obj->give_parameters());

 //        print("acceptance probs: ", 2*probn(-sqrt(T_imag/T)*sigma/sqrt(2)), 2*probn(-sqrt(T_imag/T)*kappa/2) ); 
	
  
//  DrawDensity(1,mZ ,"",  1, 1, 1);
  //	  for (decl j=0; j <10; j++)
//	  {   Draw(j, mZbig[0:i][j*10]'   ,"",  1, 1, 1);
//	   }    ShowDrawWindow();  
//	  print(  mZ[][i]);

  decl M = 10;
  for (decl i=0; i < (T+1); i++) mU[i] = ranu(dim, M+1);

  
  PF_loop_fulladapt_U(filt_obj,out_obj, g_mY, M ,  1,  mU);
  decl mExpx = out_obj->give_mExp()	;
//  Draw(0, mExpx | mY ,"",  1, 1, 1);
//   ShowDrawWindow();
//   exit(1);
 */	  
//
//  	  parameters[0][0] = transform_back(vPar)[1][];	  //sig_eta, mu, phi,
//	  parameters[1][0] = transform_back(vPar)[2][];
//  	  parameters[2][0] = transform_back(vPar)[0][];
//	  sig_eta = parameters[2][0], mu = parameters[0][0], mT = parameters[1][0]  ;
//	  mIP =  sig_eta ^2/(unit(1) - mT^2);	//	  print( sig_eta~mu~mT);
  
















	 



	 
 