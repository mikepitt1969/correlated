
% simple random effects model 
% x \sim N(theta,1),  y \sim N(x,1)
% theta_true=0.5
% Marginally,  y \sim N(theta_true,2)

clear all
T=8192;
theta_true=0.5;
y=theta_true.*ones(1,T)+sqrt(2).*randn(1,T);


%save('data','T','y','theta_true')


%save data.mat  T y theta_true

save data
