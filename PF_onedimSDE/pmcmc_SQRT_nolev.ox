// 27/10/16
// Do the one dim Nelson SV 
// This is the three parameter version - no leverage - for speed

#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include "ar1noise_class_vol.ox"
decl g_mY,  g_mSig_eps;

/*
	CPM for the SQRT process applied to no leverage case
	s&P returns until June 29, 2006.
*/

log_prior2(const parameters);
// The correlated PM (CPM) with the correlation given by phi = exp(-delta)
particle_MCMC(const aPart_output, const N , const M_eul,  const delta, const RUN_g, const BURNIN, const output_file)	 ;

main()
{
/* SHORT DATA */
	
  	decl T = 2000, SIM = 20000, dim = 1;   // state dim
	decl beta = 0.47, psi = 1  ;  //0.27, 1
  	decl N = int(beta *sqrt(T)), delta = psi * N/T ; // Euler discretisation points	(integer)
	decl  M_eul = 4, Eul_delta = 1/M_eul; 
  	decl phi = exp( -  delta), alpha ;
  	decl theta = 0.98,  mu = 1.25, sig_eta =  .142, rho = -0.0;// Nelson model parameters
  
  	
	decl mS = zeros(1, T+1);
	mS =  100 * log(loadmat("tablesp500_two.mat")'[0][(9000 - T ) :9000]) ;
	decl mY = mS[0][1:(T)] - mS[0][0:(T-1)];	// cont compounded returns
    g_mY = mY, T = columns(mY);
  
  	decl parameters = ones(4,1);  // only 3 actually estimated here
  	parameters[0][0] = mu, parameters[1][0] = theta,		  // starting values
  	parameters[2][0] = sig_eta, parameters[3][0] = rho;
	
  	decl filt_obj = new Pfilter_ar1_Data(parameters, dim, Eul_delta);  // set up PF	
  	decl out_obj = new Pfilter_ar1_output(T, dim) ;		// set up PF output
   	
	decl mU = new array[T+1]; for (decl t=0; t < (T+1); t++) mU[t] = rann(M_eul, N+1); // al the Gaussian u's
	decl mZ = zeros(1, SIM) ;	decl kappa, sigma; ranseed(100);//2.4;
	decl RUN_g = 300, BURN_IN = 3000; 	// run IID from g(U) with fixed parameters
	
	

	decl particle_output = zeros(7, SIM);
	particle_MCMC(&particle_output, N, M_eul, delta, RUN_g, BURN_IN, "output_nolev.mat");	
 	// output file can be analysed in analyse_SDE.ox
}
 
 
log_prior2(const parameters)
{
  decl mu =	parameters[0][0] , theta  =	parameters[1][0],
  sig_eta= parameters[2][0]  , rho= parameters[3][0]  ;
  return -0.5 * (mu-1)^2/0.2^2 -0.5 * (theta - 0.97)^2/0.005^2
    -0.5 * (log(sig_eta^2) - log(0.03))^2/0.1 -0.5 * (rho +0.45)^2/0.1^2;

}


 particle_MCMC(const aPart_output, const N , const M_eul,  const delta, const RUN_g, const BURN_IN, const out_file)
{
    decl phi =0, alpha;
	decl SIM = columns(aPart_output[0]);
	decl T = columns(g_mY);
	
	decl dim = 4;
//	decl vPar = zeros(dim, SIM), reject = 0.0, dof = 10;
	decl theta = 0.98,  mu = 1.25, sig_eta =  .142, rho = -0.0;// Nelson model parameters
	decl parameters = ones(4,1);  // only 3 actually estimated here
  	parameters[0][0] = mu, parameters[1][0] = theta,		  // starting values
  	parameters[2][0] = sig_eta, parameters[3][0] = rho;

	decl ratio, curr_loglik ,new_loglik, sigma, kappa ;
	decl logMH=0, parameters_new = parameters, logL = -10e6, logL_new, accept=0,
	log_prnew, log_pr =0;

	decl mP = zeros(4, SIM);
	decl mZ = zeros(2, SIM);	  // decl  = zeros(7, SIM);

	decl mU = new array[T+1]; for (decl t=0; t < (T+1); t++) mU[t] = rann(M_eul, N+1); // al the Gaussian u's
	decl mU_new = mU;
		decl filt_obj = new Pfilter_ar1_Data(parameters, dim, 1/M_eul);  // set up PF	
  	decl out_obj = new Pfilter_ar1_output(T, dim) ;		// set up PF output
	decl mY = g_mY;
	for (decl j=0; j<SIM; j++)
	{
	
        parameters_new[0][0] = exp(log(parameters[0][0]) + 3/sqrt(T) * rann(1,1));	 // mu
		parameters_new[1][0] =  0.996 * probn( quann(parameters[1][0]/0.996)+ 3/sqrt(T) * rann(1,1));	 // theta - pers
		parameters_new[2][0] = exp(log(parameters[2][0]) + 2.5/sqrt(T) * rann(1,1) ) ;   //2/sqrt(T), sigma
		parameters_new[3][0] = 0; 	// leverage fixed at zero for now
						 // 2 * probn( quann(0.5 * parameters[3][0]+0.5)+ 3/sqrt(T) * rann(1,1) ) -1;

		if ( j <= BURN_IN) parameters_new = parameters	;
		alpha = 2 * parameters_new[0][0] * (1-	parameters_new[1][0])/ parameters_new[2][0]^2;	// marginal dof in gamma process
		if (alpha < 1) logMH = -1000;
		else
		{
			filt_obj->reset_parameters(parameters_new);
	  		for (decl t=0; t < (T+1); t++)
	  		{
	  			mU_new[t] = phi .* mU[t] + sqrt(1-phi.^2) .* rann(M_eul,N+1);   //	probn( )
	  		}
			PF_loop_U_sort(filt_obj, out_obj, mY,	N,  1, mU_new)	;
			logL_new = sumr( out_obj->give_mlogf() ) ;
			log_prnew = log_prior2(parameters_new);  
			logMH  = logL_new - logL + log_prnew - log_pr;
		}
		phi = exp(-delta);
		if (j <= RUN_g) {logMH = 0; phi=0;}
		mZ[0][j] = logL_new - logL; // W really when param fixed
		if (isnan(logMH)) {
			print( logL	~ logL_new~alpha , parameters~parameters_new );
			print(out_obj->give_mlogf()[][]' );
			exit(1);
		}
		if (log(ranu(1,1)) < logMH)		 //	 logMH
		{
		   for (decl t=0; t < (T+1); t++) mU[t] = mU_new[t];
		   logL = logL_new;
		   log_pr = log_prnew;
		   parameters = parameters_new ;
		   accept+=1;
		}			   
		
	  	mZ[1][j] = logL;	mP[][j]	= parameters[][0] ;
//		particle_output[][j] = 
		aPart_output[0][][j] = mZ[0][j] | logL	| accept/(j+1) | parameters[][0] ;
			if (j > RUN_g) 	aPart_output[0][2][j] = (accept-RUN_g) /(j-RUN_g+1);
		if (j > RUN_g && fmod(j, 500) ==0 )
		{
			sigma = sqrt( varr( mZ[1][0:(RUN_g-1)] ) );
		    kappa =  sqrt( varr(mZ[0][RUN_g+2:j]) ) ;
			print("i, N, T",j~N~T~logL~logMH~accept/j);
			print("sigma2", sigma^2);
			print("mean, SD Z, proba, delta", meanr(mZ[0][100:j])~kappa~2 * probn(-kappa/2)~delta  );
			kappa = sqrt(-2 *meanr(mZ[0][1:j]) );
			print("kappa mean", 2*probn(-kappa/2) ) ; 
			Draw(0, mZ[0][1:j]  ,"",  1, 1, 1);
			Draw(1, mZ[1][0:j]  ,"",  1, 1, 1);
			for (decl jj =0; jj < 4; jj++) Draw(2+jj, mP[jj][0:j]  ,"",  1, 1, 1);
			savemat(out_file,aPart_output[0]);
	    	ShowDrawWindow();
		}
	 }
	 savemat(out_file,aPart_output[0]);
}








	 



	 
 