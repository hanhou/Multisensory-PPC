
function [info_out info_train_out w_ole est_angle perc_right_test] =...
          compute_info(data1,data2,step_test)
global Y_targ X_train

aux = size(data1);
nu_trial = aux(2)/2;
N_out = aux(1);

Y_targ(1,1:nu_trial) = 1;
Y_targ(1,nu_trial+1:nu_trial*2) = -1;


X_train = data1 - mean(data1')'*ones(1,nu_trial*2);
X_test  = data2 - mean(data2')'*ones(1,nu_trial*2);
X_train(N_out+1,:)= ones(1,nu_trial*2);
X_test(N_out+1,:)= ones(1,nu_trial*2);

[info_out info_train_out w_ole] = cg(X_test,N_out,nu_trial,step_test);
                                  
est_angle = w_ole'*X_test;
perc_right_test = sum(max(Y_targ)*((est_angle>0)*2-1)==Y_targ)/(nu_trial*2)*100;


