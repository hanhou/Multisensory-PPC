clear 
clear global

load data_out2
load data_in2
load w_ole_out
load w_ole_in

step_test = 1;

aux = size(data_out);
N_out= aux(1);
nu_trial = aux(2)/2;


aux = size(data_in);
N_in= aux(1);

Y_targ(1,1:nu_trial) = 1;
Y_targ(1,nu_trial+1:nu_trial*2) = -1;

data_out(N_out+1,:)= ones(1,nu_trial*2);
est_angle_out = w_ole_out'*data_out;
perc_right_test_out = sum(max(Y_targ)*((est_angle_out>0)*2-1)==Y_targ)/(nu_trial*2)*100;
info_out = fisher(w_ole_out,data_out',step_test);

data_in(N_out+1,:)= ones(1,nu_trial*2);
est_angle_in = w_ole_in'*data_in;
perc_right_test_in= sum(max(Y_targ)*((est_angle_in>0)*2-1)==Y_targ)/(nu_trial*2)*100;
info_in = fisher(w_ole_in,data_in',step_test);

fprintf('Number of trial        : %d\n',nu_trial);

fprintf('Info In test           : %2.3f\n',info_in);
fprintf('Percentage classification correct: %2.2f\n',perc_right_test_in);
fprintf('\n');

fprintf('Info Out test          : %2.3f\n',info_out);
fprintf('Percentage classification correct: %2.2f\n',perc_right_test_out);
fprintf('\n');