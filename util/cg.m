%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information estimated with conjugate gradients methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [info_test info_train w_ole_cg] = cg(X_test,nu_unit,nu_trial,step_test);
global Y_targ X_train

fprintf('Initializing weights\n');
w_ole_ini = rand(nu_unit+1,1)*1e-8; 
w_ole_ini = w_ole_ini-mean(w_ole_ini);
w_ole_ini(nu_unit+1)=0; 


disp('Minimization using conjugate gradients methods');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X_train must be Nu_unit x (2*Nu_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux= size(X_train);

if (aux(1)==nu_unit+1) & (aux(2)==2*nu_trial)
    disp('Input data has the right format');
else
    error('Error in cg.m: input data has the wrong format')
end


gtol = eps; 
ftol=1e-14;
iters=1;
c1=0.001; c2=0.6; %%necessary 0<c1<c2<1 Wolfe line search
info_test_min=0;
count=0;
count2=0;
k=0;
max_iter = 10000;
last_info = 1e20;

while k<max_iter
	   [w_ole_cg,err,foo1,foo2,foo3,exitflag] = ...
           cg_pr('foncgradOLE2','gradOLE2',w_ole_ini,gtol,ftol,c1,c2,iters,0); % Polak-Ribiere conjugate gradient method (without restart) 
 info_train = fisher(w_ole_cg,X_train',step_test);
 info_test = fisher(w_ole_cg,X_test',step_test);
 fprintf('Iter: %d  I_train: %6.6f  I_test: %6.6f  I_test max: %6.6f  c:%d  c2:%d\n',...
          k*iters,info_train,info_test,info_test_min,count,count2);      

   if abs(info_test-last_info)<last_info*0.001
      count2 = count2+1;
      if count2==3
         k=max_iter;
      end
   end

   if info_test_min<info_test
      info_test_min=info_test;
      w_ole=w_ole_cg;
      count=0;
   else
      count=count+1;
      if count==20
         k=max_iter;
      end
   end    
   k=k+1;
   w_ole_ini=w_ole_cg;
   last_info = info_test;
end


theta_hat =  w_ole_cg'* X_train;
mean_theta_hat1 = mean(theta_hat(1:nu_trial));
mean_theta_hat2 = mean(theta_hat(nu_trial+1:2*nu_trial));
var_theta_hat1 = var(theta_hat(1:nu_trial));
var_theta_hat2 = var(theta_hat(nu_trial+1:2*nu_trial));
slope = (mean_theta_hat1-mean_theta_hat2)/step_test;
std_av = ((var_theta_hat1+var_theta_hat2)/2)^0.5;
information = (slope/std_av)^2;

%%% properly normalize the weights
%%% w_ole_cg at this point is h' 
w_ole_cg=w_ole_cg*abs(slope)/std_av^2;

info_train = fisher(w_ole_cg,X_train',step_test);
