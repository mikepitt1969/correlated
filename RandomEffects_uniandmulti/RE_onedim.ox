#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <quadpack.h>
#pragma link("maximize.oxo")

decl g_mY;
decl g_sigma;
/* Function declarations - functions defined at the bottom */
// For IACT estimation
Parzen( const lag);
// mX 1*T.
// IACT estimation
var_inflate(const lag, const mX);  
// The simulated estimator for the log-lik
est_lik(const vPar,const mX_tilda,  const N);
// derivate - not used here
est_weights_dash(const vPar,const mX_tilda,  const N) ;
// the true likelihood - as it is a Gaussian RE model
tru_lik(const vPar);
// Used for maximisation - when necessary
lik_max(const vPar, const adFunc, const avScore, const amHessian);
// Exact MCMC - in other words using the known true likelihood and a simple RW (t-distribution) proposal 
exact_MCMC(const aExact_output, const mu_mean, const mu_var);
// The correlated PM (CPM) with the correlation given by rho
// Same proposal (in theta) as for the Exact MCMC routine above
particle_MCMC(const aPart_output, const mu_mean, const mu_var,
		 const N, const rho);

/* main function */

main() {
   decl T = 8192/8, SIM = 32000;
   decl beta = 0.7, psi = 0.7; // for sim 1
   decl N = int(beta * T^0.5), rho = exp(-psi * N/T) ;	
   
// Model and model generation:
	ranseed(12);	  //14, 12   
   decl  mu = 0.5, sigma = sqrt(1); g_sigma = sigma; // sqrt(1)
   decl mX = mu +  rann(1, T);
   decl mY = mX + sigma .* rann(1, T); g_mY = mY;
// Get the curvature for RW proposal var:
   decl vPar = zeros(1,1); vPar[0][] = mu;  
   decl mu_mean = 0.5, mu_var =   unit(2);
   Num2Derivative(lik_max, mu_mean, &mu_var);  mu_var = - 1/(mu_var);  //use for RW proposal var

   decl lag_EX = 60, lag_PF = 220, Brn = 2000 ;
   decl mExact_output = zeros(1, SIM), mPart_output = zeros(4, SIM)	;
/* EXACT MCMC using the true likelihood and RW proposal */
   decl EX_ACC = exact_MCMC(&mExact_output, mu_mean, mu_var);
   decl EX_IF  = var_inflate(lag_EX, mExact_output[0][Brn:(SIM-1)]) ;

   ranseed(90);
/* CPM using the estimated likelihood and RW proposal */
   decl PF_ACC = particle_MCMC(&mPart_output, mu_mean, mu_var,  N,  rho)	;
   vPar = 	mPart_output[0][0:(SIM-1)];	// theta draws
   decl mZp =  mPart_output[1][0:(SIM-1)];	// proposed Z
   decl mZa =  mPart_output[2][0:(SIM-1)] ;	//accepted Z

/* Display the output */
   Draw(0,  mZa,0, 1);DrawTitle(0, " Accepted Z" );	   // accepted Z's  
   Draw(1,  mZp-mZa,0, 1); DrawTitle(1, " R=W - Z" ); // difference Z'-Z
   mZa = mZa[][Brn:(SIM-1)], mZp = mZp[][Brn:(SIM-1)] ;
   decl s2 = varr(mZa);
   DrawDensity(2, mZa  ," Accepted Z marginal",  1, 1, 1);	//accepted
   DrawXMatrix(2,  1/sqrt(s2).*densn((sortr(mZa) -s2/2)/sqrt(s2)), {" Theoretical approx "}, sortr(mZa) , "", 0, 1);
   decl kappa2 = varr(mZp - mZa);
   DrawDensity(3, mZp - mZa ," R=W - Z",  1, 1, 1); // R = Z' -Z
   DrawXMatrix(3,  1/sqrt(kappa2).*densn((sortr(mZp - mZa) +kappa2/2)/sqrt(kappa2)), {" Theoretical approx "}, sortr(mZp - mZa) , "", 0, 1);	
   Draw(4, vPar,0, 1);	DrawTitle(4, " theta " ); 
   DrawCorrelogram(5,vPar[0][Brn:(SIM-1)] , "Correlogram  theta  ", lag_PF);
   DrawCorrelogram(6, mZp - mZa, " Correlogram R = Z'-Z", lag_PF);	   //correlogram Z' - Z
   //SaveDrawWindow("out.eps");
   ShowDrawWindow();
/* Write a few results to console */
   
   decl Pr_th =  2 * probn(-sqrt(kappa2)/2)	;	 // theor prob
   decl PF_IF = var_inflate(lag_PF, mPart_output[0][Brn:(SIM-1)]);  	 

   print("Summary: rho, N, T: ", rho~N~T, "\n");
   print("Summary: Exact acceptance and ineff, CPM acceptance and inefficiency: ", EX_ACC~EX_IF~PF_ACC~PF_IF, "\n");
   print("mean and var of Z (accepted values)", meanr(mZa)~varr(mZa ));
   print("mean and var of R=W-Z", meanr(mZp - mZa)~varr(mZp - mZa ) );
   print("Theor Lower Bnd on acceptance prob and actual acceptance prob: ", Pr_th * EX_ACC ~ PF_ACC);    

}


// For IACT estimation
Parzen( const lag)
{
  decl working, i;
  decl parzen = zeros(1, lag);
  for (i = 0; i < lag; i++)
  {
      working = double(i)/double(lag);
      if (working < 0.5)		  
		parzen[][i] = 1.0 - (6.0 * (working)^2) + (6.0 * (working)^3);
      else
		parzen[][i] = 2.0 * (1.0 - working)^3 ;
  }
  return parzen;
}

// mX 1*T.
// IACT estimation
var_inflate(const lag, const mX)
{ 
  decl acfX = acf(mX', lag)[1:lag][];
 
  decl par = Parzen(lag);
  return 1 + 2 * sumc(acfX ); //.* par'
}

  
// The simulated estimator for the log-lik
est_lik(const vPar,const mX_tilda,  const N)
{
   decl mY = g_mY; // 1 * T
   decl T = columns(mY);
   decl mu = vPar[0][] ;
   decl sigma= g_sigma; // meas noise - global variable
   decl mX = mu + mX_tilda;	   // N * T
   decl mA = densn((mY - mX)./sigma)./(N *sigma);	  // N * T
   decl mAv = sumc(mA);	  //sumr over N	 - get 1*T vector of averages
   decl logL = log(mAv);  // 1 * T 
   return sumr( logL );   //sum over T for log-lik estimate
}
// derivate - not used here
est_weights_dash(const vPar,const mX_tilda,  const N)
{
   decl mY = g_mY; // 1 * T
   decl T = columns(mY);
   decl mu = vPar[0][] ;
   decl sigma= g_sigma; //= exp(vPar[1][])  ;   
   decl mX = mu + mX_tilda;
   decl mA1 = densn( (mY - mX)./sigma ) .* ((mY - mX)./sigma) ./(N *sigma^2);
   decl mAv1 =  sumc(mA1.^2);
   return ( mAv1); 
}

// the true likelihood - as it is a Gaussian RE model
tru_lik(const vPar)
{
   decl mY = g_mY; // 1 * T
   decl mu = vPar[0][] ;
   decl sigma = g_sigma; //exp(vPar[1][])  ;
   decl logL = -0.5 * log(2*3.142) -0.5 * log(1+sigma^2)
        -0.5 * ((mY - mu).^2) ./(1+sigma^2);
   return sumr( logL );	  // sum over T
}
// Used for maximisation - when necessary
lik_max(const vPar, const adFunc, const avScore, const amHessian)
{  	adFunc[0] = (tru_lik(vPar[][]))-0.5 * ((vPar[0][] - 0.5)^2)/0.25 ;  return 1;} 


// Exact MCMC - in other words using the known true likelihood and a simple RW (t-distribution) proposal 
exact_MCMC(const aExact_output, const mu_mean, const mu_var) {
	
	decl SIM = columns(aExact_output[0]);
	decl ratio, curr_loglik = -10e6 ,new_loglik ;
	decl vPar = zeros(1, SIM), reject = 0.0, dof = 10;
	decl pr_mean = 0.5, pr_var = 0.25;
	decl new_prior, curr_prior; 
	vPar[0][0] = 0.5;
	for (decl i=1 ; i < SIM; i++)
	{  
		vPar[][i] = vPar[][i-1] +  sqrt(mu_var)* rant(1,1, dof);
		new_prior = -0.5 * ((vPar[][i] - pr_mean)^2)/pr_var;
		curr_prior  = -0.5 * ((vPar[][i-1] - pr_mean)^2)/pr_var; 
		new_loglik = tru_lik(vPar[][i])	;	 
		ratio = new_loglik + new_prior - curr_loglik - curr_prior	; // log-lik terms and prior
		
		if (ranu(1,1) > exp(ratio)) {
			vPar[][i] =	vPar[][i-1] ;	reject += 1; //rejected
		}
		else  curr_loglik = new_loglik;
	}
	aExact_output[0] = 	vPar ;
	print("acceptance rate standard MCMC: ", 1-reject/SIM, "\n"); // av acceptance
	return 1-reject/SIM;
}

// The correlated PM (CPM) with the correlation given by rho
// Same proposal (in theta) as for the Exact MCMC routine above
particle_MCMC(const aPart_output, const mu_mean, const mu_var,
		 const N, const rho)
{	
	decl SIM = columns(aPart_output[0]);
	decl T = columns(g_mY);
	decl ratio, curr_loglik ,new_loglik ;
	decl vPar = zeros(1, SIM), reject = 0.0, dof = 10;
	decl mX_tilda_curr = zeros(N,T)	;
	decl mX_tilda_new = zeros(N,T)	; vPar[0][0] = 0.5;
	mX_tilda_curr = rann(N,T);
	decl mZ_prop = zeros(1, SIM), mZ_acc = zeros(1, SIM);
	decl pr_mean = 0.5, pr_var = 0.25;
	decl new_prior, curr_prior; 
	curr_loglik = est_lik(vPar[][0], mX_tilda_curr, N) //tru_lik(vPar[][i])
			-0.5 * ((vPar[][0] - pr_mean)^2)/pr_var	;
			
	for (decl i=1 ; i < SIM; i++)
	{	 
		vPar[][i] = vPar[][i-1] + sqrt(mu_var)* rant(1,1, dof);  // RW proposal	 	
		// The correlated proposal for the Gaussian terms
		mX_tilda_new = rho .* mX_tilda_curr + sqrt(1-rho^2).* rann(N,T);
		// Done 
		new_loglik =  est_lik(vPar[][i], mX_tilda_new, N);
		new_prior = -0.5 * ((vPar[][i] - pr_mean)^2)/pr_var; //with prior
		curr_prior = -0.5 * ((vPar[][i-1] - pr_mean)^2)/pr_var	 ;

		ratio = new_loglik + new_prior - curr_loglik - curr_prior;

		mZ_prop[][i] = new_loglik - tru_lik(vPar[][i]);  // Not need for alg - just to record proposed Z
		mZ_acc[][i] = curr_loglik - tru_lik(vPar[][i-1]);	// Just record accepted Z's 

	
		if (ranu(1,1) > exp(ratio)) {
			vPar[][i] =	vPar[][i-1] ; reject += 1; //rejected
		}
		else  {
			curr_loglik = new_loglik;	mX_tilda_curr = mX_tilda_new;  //accepted
		}
	    if (fmod(i, 500) ==0 ) {
				print("acc rate PF: ",i~ 1-reject/(i+1));
				
		}
	
	}
	aPart_output[0][0][] = vPar ;
	aPart_output[0][1][] = mZ_prop;
	aPart_output[0][2][] = mZ_acc;
	print("acceptance rate PMCMC: ",  1-reject/SIM, "\n");
	return 1-reject/SIM;
}



