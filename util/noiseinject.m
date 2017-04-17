clear
clear global

global Y_targ X_train;
rand('state',sum(100*clock))

N_in = 200;
pos_in=[0:180/N_in:179.9999];
N_out = 200;
pos_out=[0:180/N_out:179.9999];
 
peak = 90;
trial_dur = 500;
delta_t = 1e-3;
 
step_test =3;
nu_trial = 10001;
max_rate_in = 10;
K_in = 2;
K_oin = 5;
K_oo =20;
K_oI = 3;
threshold = 0.0;
 
length_ker = 0;
time_const_out = 10e-3;
time_const_in = 10e-3;
 
g_w_oi = 9.0;
g_w_oo = 7.0;
amp_i = 0.9; 



%amplitude of inhibitory surround for w_oo

shuffle_flag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes_in = zeros(N_in,trial_dur+length_ker);
spikes_in_trial = zeros(N_in,nu_trial*2);
spikes_in2 = zeros(N_in,trial_dur+length_ker);
spikes_in_trial2 = zeros(N_in,nu_trial*2);

aux1_proba1 = zeros(N_out,trial_dur+2);
aux2_proba1 = zeros(N_out,trial_dur+2);
aux1_proba2 = zeros(N_out,trial_dur+2);
aux2_proba2 = zeros(N_out,trial_dur+2);

proba_out  = zeros(N_out,trial_dur+length_ker);
spikes_out = zeros(N_out,trial_dur+length_ker);
spikes_out_trial  = zeros(N_out,nu_trial*2);
proba_out_trial  = zeros(N_out,nu_trial*2);
proba_out2  = zeros(N_out,trial_dur+length_ker);
spikes_out2 = zeros(N_out,trial_dur+length_ker);
spikes_out_trial2 = zeros(N_out,nu_trial*2);
proba_out_trial2  = zeros(N_out,nu_trial*2);

w_oi = zeros(N_out,N_in);
w_oo = zeros(N_out,N_out);

for j=1:N_out
    w_oi(j,:) = g_w_oi/N_in *(exp(K_oin*(cos((pos_in-pos_out(j))/180 *2*pi)-1)));
    w_oo(j,:) = g_w_oo/N_out*...
                ((exp(K_oo*(cos((pos_out-pos_out(j))/180*2*pi)-1)))-...
           amp_i*(exp(K_oI*(cos((pos_out-pos_out(j))/180*2*pi)-1))));
end


ker_in = exp(-([length_ker:-1:1]-1)/time_const_in)';
ker_out = exp(-([length_ker:-1:1]-1)/time_const_out)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%
count =1;
for l=1:2
   pos = peak+step_test*(l-1.5);
   proba_in = max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))*delta_t; 
   fprintf('%d  pos: %4.2f\n  ',l, pos);
   for j=1:nu_trial   
       if round(j/50)==j/50
          fprintf('%d  ',j);
       end

      spikes_in(:,1:trial_dur)  = rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));
      spikes_in2(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));
                   
      for k=1:trial_dur-1
         aux1_proba1(:,k+1) = (1-delta_t/time_const_in)*aux1_proba1(:,k)...
                              +1/time_const_in*(w_oi*spikes_in(:,k));
         aux2_proba1(:,k+1) = (1-delta_t/time_const_out)*aux2_proba1(:,k)+...
                              +1/time_const_in*(w_oo*spikes_out(:,k));
                            
         proba_out(:,k+1) = aux1_proba1(:,k+1) + aux2_proba1(:,k+1);
                            
         aux1_proba2(:,k+1) = (1-delta_t/time_const_in)*aux1_proba2(:,k)...
                              +1/time_const_in*(w_oi*spikes_in2(:,k));
         aux2_proba2(:,k+1) = (1-delta_t/time_const_out)*aux2_proba2(:,k)...
                              +1/time_const_in*(w_oo*spikes_out2(:,k));
         proba_out2(:,k+1) = aux1_proba2(:,k+1) + aux2_proba2(:,k+1);
                    
         spikes_out(:,k+1)  = (rand(N_out,1) < delta_t*(proba_out(:,k+1) -threshold));
         spikes_out2(:,k+1) = (rand(N_out,1) < delta_t*(proba_out2(:,k+1)-threshold));
      end
   
      proba_out_trial(:,count)=sum(proba_out');
      proba_out_trial2(:,count)=sum(proba_out2');
      
      spikes_in_trial(:,count) = sum(spikes_in')';
      spikes_out_trial(:,count) = sum(spikes_out')';
      spikes_in_trial2(:,count) = sum(spikes_in2')';
      spikes_out_trial2(:,count) = sum(spikes_out2')';
      count = count+1;
   end
   fprintf('\n');   
end


%%%%%%%%%%%%%%%%%%%%%%%%
% stats
%%%%%%%%%%%%%%%%%%%%%%%%
mean_in = mean(spikes_in_trial(:,1:nu_trial)')/(trial_dur*delta_t);
mean_in2 = mean(spikes_in_trial(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);
var_in = var(spikes_in_trial(:,1:nu_trial)')/(trial_dur*delta_t);

mean_out = mean(spikes_out_trial(:,1:nu_trial)')/(trial_dur*delta_t);
mean_out2 = mean(spikes_out_trial(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);

mean_proba  = mean(proba_out_trial(:,1:nu_trial)');
mean_proba2 = mean(proba_out_trial(:,nu_trial+1:2*nu_trial)');

var_out = var(spikes_out_trial(:,1:nu_trial)')/(trial_dur*delta_t);
cov_in =  cov(spikes_in_trial(:,1:nu_trial)')/(trial_dur*delta_t);
cov_out =  cov(spikes_out_trial(:,1:nu_trial)')/(trial_dur*delta_t);

difference = mean(spikes_out_trial(:,1:nu_trial)')...
            -mean(spikes_out_trial(:,nu_trial+1:nu_trial*2)');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rate info
Y_targ(1,1:nu_trial) = 1;
Y_targ(1,nu_trial+1:nu_trial*2) = -1;

X_train = proba_out_trial - ...
                    mean(proba_out_trial')'*ones(1,nu_trial*2);
X_test = proba_out_trial2 - ...
                    mean(proba_out_trial2')'*ones(1,nu_trial*2);
X_train(N_out+1,:)= ones(1,nu_trial*2);
X_test(N_out+1,:)= ones(1,nu_trial*2);

[info_rate_out info_rate_train_out w_ole_rate_out] = cg(X_test,N_out,nu_trial,step_test);

est_rate_angle = w_ole_rate_out'*X_test;
perc_right_rate_test= sum(max(Y_targ)*((est_rate_angle>0)*2-1)==Y_targ)/(nu_trial*2)*100;
    
    


% output info
Y_targ(1,1:nu_trial) = 1;
Y_targ(1,nu_trial+1:nu_trial*2) = -1;

spikes_out_trial = spikes_out_trial - ...
                    mean(spikes_out_trial')'*ones(1,nu_trial*2);
spikes_out_trial2 = spikes_out_trial2 - ...
                    mean(spikes_out_trial2')'*ones(1,nu_trial*2);
X_train = spikes_out_trial;
X_train(N_out+1,:)= ones(1,nu_trial*2);
X_test = spikes_out_trial2;
X_test(N_out+1,:)= ones(1,nu_trial*2);

[info_out info_train_out w_ole_out] = cg(X_test,N_out,nu_trial,step_test);

est_angle = w_ole_out'*X_test;
perc_right_test= sum(max(Y_targ)*((est_angle>0)*2-1)==Y_targ)/(nu_trial*2)*100;
    
    
%input info
spikes_in_trial = spikes_in_trial - ...
                    mean(spikes_in_trial')'*ones(1,nu_trial*2);
spikes_in_trial2 = spikes_in_trial2 - ...
                    mean(spikes_in_trial2')'*ones(1,nu_trial*2);
                
X_train = spikes_in_trial;
X_train(N_in+1,:)= ones(1,nu_trial*2);
X_test = spikes_in_trial2;
X_test(N_in+1,:)= ones(1,nu_trial*2);
[info_in info_train_in w_ole_in] = cg(X_test,N_in,nu_trial,step_test);

est_angle_in = w_ole_in'*X_test;
perc_right_test_in = sum(max(Y_targ)*((est_angle_in>0)*2-1)==Y_targ)/(nu_trial*2)*100;
 
%normalize for duration of trial
info_in = info_in/(trial_dur*delta_t);
info_out= info_out/(trial_dur*delta_t);
info_rate_out= info_rate_out/(trial_dur*delta_t);
info_train_in = info_train_in/(trial_dur*delta_t);
info_train_out= info_train_out/(trial_dur*delta_t);
info_rate_train_out= info_rate_train_out/(trial_dur*delta_t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute true info 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:10
    
   pos = peak;
   proba_in = max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))*delta_t; 
   
   noise_image = j;

   ac_in = max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1));
   ac_in_der = -max_rate_in * K_in * sin((pos_in'-pos)/180*2*pi).*...
               exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))/180*2*pi;

   info_in_true = sum(ac_in_der.^2./(ac_in+noise_image));



true_cov_in =diag(ac_in);

local_cov_out = diag(mean_out2);

% 1e-6 is added to ensure invertibility in the next line when aux_info is
% zero
aux_info = w_oo*local_cov_out*w_oo' + diag(diag(ones(N_out)*1e-6));
info_out_rate_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);

aux_info = local_cov_out + diag(diag(ones(N_out)*1e-6));
info_out_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);

end

%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n');
fprintf('True Info In test       : %2.3f\n',info_in_true);
fprintf('True Info rate Out test : %2.3f\n',info_out_rate_true);
fprintf('True Info Out test      : %2.3f\n',info_out_true);
fprintf('\n');


fprintf('Info rate Out test          : %2.3f\n',info_rate_out);
fprintf('Info rate Out train         : %2.3f\n',info_rate_train_out);
fprintf('Percentage classification correct: %2.2f\n',perc_right_rate_test);
fprintf('\n');

fprintf('Max firing rate: %2.3f\n',max(mean_out));
fprintf('\n');


subplot(221)

subplot(222)


subplot(223)



subplot(224)
