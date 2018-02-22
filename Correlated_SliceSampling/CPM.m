% Correlated pseudo marginal
% Basic random effect models X \sim N(theta,1),  Y \sim N*(X,1) with
% theta_true=0.5


load data                   % data are simulated using simul_data.m 

niter=10000;                % number MCMC iterations
nburnin=2000;               % number of burnin iterations

N=56;                       % number of particles
beta=0;                     % parameter of correlation coef for proposal U
psi=0;                      % parameter of correlation coef for proposal U 
pho=0.9962;                 % correlation coefficient for proposal U: pho=exp(-psi*beta/sqrt(T))
std_prop=1./sqrt(T);        % standard proposal distribution theta


thetas=zeros(1,niter);
llikes=zeros(1,niter);

llike_cand=0;               % candidate simulated log-likelihood
u_cand=zeros(1,T*N);        % candidate auxiliary variate


theta_current=theta_true;
u_current=randn(1,T*N);
llike_current=llikelihood(y,theta_current,u_current);

thetas(1,1)=theta_current;
llikes(1,1)=llike_current;

tic 

for i=2:niter
    
    u_cand=pho.*u_current+sqrt(1-pho^2).*randn(1,T*N);
    
    theta_cand=theta_current+std_prop*randn(1);
    
    llike_cand=llikelihood(y,theta_cand,u_cand);
     
    if log(rand(1))<(llike_cand-llike_current)       % random walk proposal + flat prior on theta
        
        u_current=u_cand;
        
        llike_current=llike_cand;
        
        theta_current=theta_cand;
        
    end
        
thetas(1,i)=theta_current;
llikes(1,i)=llike_current;
        
end

acf_theta=autocorr(thetas(1,nburnin:end),50);
acf_llike=autocorr(llikes(1,nburnin:end),50);

plot(acf_theta);

save results_CPM_T=8192 thetas llikes T y N acf_theta acf_llike

toc
    
