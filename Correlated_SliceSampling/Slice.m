% Pseudo-marginal Slice Sampling
% Update parameter theta using RW proposal given auxiliary variates U
% Update auxiliary variates U given theta using Elliptical Slice Sampling

% Basic random effect models X \sim N(theta,1),  Y \sim N*(X,1) with
% theta_true=0.5


load data                   % data are simulated using simul_data.m 

niter=10000;                % number MCMC iterations
nburnin=2000;               % number of burnin iterations

N=56;                       % number of particles

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

u_ess=zeros(1,T*N);         % allows to define ellipse in elliptical slice sampling
angle_ess=0;
llike_ess=0;
angle_min=0;
angle_max=0;

tic

for i=2:niter

    if mod(i,1000)==0 
        i 
    end
    
    % update theta given u
    
     theta_cand=theta_current+std_prop*randn(1);
    
     llike_cand=llikelihood(y,theta_cand,u_current);

     if log(rand(1))<(llike_cand-llike_current)
         
        llike_current=llike_cand;
        
        theta_current=theta_cand;
        
     end
     
     % update u given theta using Elliptical Slice sampling
   
     u_ess=randn(1,T*N);               
    
     llike_ess=llike_current+log(rand(1));
     
     angle_ess=2*pi*rand(1);
     
     angle_min=angle_ess-2*pi;
     
     angle_max=angle_ess;
     
     u_cand=u_current.*cos(angle_ess)+u_ess.*sin(angle_ess);
     
     llike_cand=llikelihood(y,theta_current,u_cand);
     
     while (llike_cand<=llike_ess)
         
        if (angle_ess<0)
            
            angle_min=angle_ess;
            
        else
            
            angle_max=angle_ess;
            
        end
        
        angle_ess=angle_min+(angle_max-angle_min)*rand(1);
        
        u_cand=u_current.*cos(angle_ess)+u_ess.*sin(angle_ess);
        
        llike_cand=llikelihood(y,theta_current,u_cand);    
       
     end
     
     u_current=u_cand;   % elliptical slice sampling always accept!
     llike_current=llike_cand;
     
thetas(1,i)=theta_current;
llikes(1,i)=llike_current;
        
end
    
acf_theta=autocorr(thetas(1,nburnin:end),50);
acf_llike=autocorr(llikes(1,nburnin:end),50);

plot(acf_theta);

save results_slice_T=8192 thetas llikes T y N acf_theta acf_llike

toc

