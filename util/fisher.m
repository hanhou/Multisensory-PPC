function information=fisher(w,act,step_test)

[nu_trial aux] = size(act);
nu_trial=nu_trial/2;

theta_hat =  w'* act';
mean_theta_hat1 = mean(theta_hat(1:nu_trial));
mean_theta_hat2 = mean(theta_hat(nu_trial+1:2*nu_trial));
var_theta_hat1 = var(theta_hat(1:nu_trial));
var_theta_hat2 = var(theta_hat(nu_trial+1:2*nu_trial));

slope = (mean_theta_hat1-mean_theta_hat2)/step_test;
std_av = ((var_theta_hat1+var_theta_hat2)/2)^0.5;
information = (slope/std_av)^2;


