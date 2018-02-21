% Return simulated log likelihood for parameter theta, auxiliary variates
% u, data y

function llike = llikelihood(y,theta,u)

T=length(y);
N=length(u)/T;

weight=zeros(T,N);

llike=0;

%
% for t=1:T
%     
%     for n=1:N
%         
%         weight(t,n)=exp(-(y(1,t)-theta-u(1,(t-1)*N+n))^2/2)/sqrt(2*pi);
%         
%     end
%     
%     llike=llike+log(mean(weight(t,:)));
%      
% end


for t=1:T
   
   weight(t,:)=exp(-((y(1,t)-theta).*ones(1,N)-u(1,((t-1)*N+1):(t*N))).^2./2);
   
   llike=llike+log(mean(weight(t,:))./sqrt(2*pi));
   
end
