%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
%% this version only adds the spike from MT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear
clear global

load h

global Y_targ X_train;
rand('state',sum(100*clock))

N_in = 100;
pos_in=[0:180/N_in:179.9999];
N_out = 100;
pos_out=[0:180/N_out:179.9999];
nu_targets = 8;
decis_thres = 28;

peak = 90;
pre_trial_dur= 50; %% must be 50 or more
trial_dur = 401;

if pre_trial_dur>trial_dur
   error('pre_trial_dur cannot be longer than trial_dur'); 
end

delta_t = 1e-3;
% difference in s between the two test stimuli
step_test =3;

nu_trial = 100;
coherence = 15;

% parameters for MT from Mazurek and Shadlen, except width. 
% past.
r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;

% parameter of the visually evoked activity
max_rate_in2 = 30/1.5;
b_in2 = 0;
K_in2 = 4;

% Weight parameters
b_out = 0;
K_oin = 5;
K_oin2= 5;
K_oo = 10;
K_oI = 2;
threshold = 0.0;

%time constant pre trial
pre_time_const_out = 1000e-3;
pre_time_const_in  = pre_time_const_out;

%time constant for integration
time_const_out = 1000e-3;
time_const_in  = time_const_out;
 
% MT weight 
g_w_oi = 14;
% Visual evoked weight
g_w_oi2= 4;
%%% drop in attention to visual target once motion stimulus appears.
att_gain2 = 0.1*0;


g_w_oo = 20;
amp_i = 0.0; 
dc_w_oo = -0.078;
 
ext_noise=0;

shuffle_flag = 1;
info_flag = 1;
shuffle_flag = 1;

%%% duration of the window used to compute Fisher
wind_size = 100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
count =1;
%% there are two main loops, for two different positions of the input
for l=1:2
   pos = peak+step_test*(l-1.5);
   fprintf('%d  pos: %4.2f\n  ',l, pos);
   
   %%% proba_in: proba of firing of input neurons in response to motion
   proba_in = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
   
   %%% proba_in2: proba of firing of input neurons in response to visual
   %%% targets
   pos_targ = pos+[0:180/nu_targets:179];
   proba_in2 = zeros(N_out,1);
   for m=1:nu_targets
       aux(:,m)= ((max_rate_in2-b_in2)*exp(K_in2*(cos((pos_in'-pos_targ(m))/180*2*pi)-1))+b_in2)*delta_t;
       proba_in2 = proba_in2 + ... 
       ((max_rate_in2-b_in2)*exp(K_in2*(cos((pos_in'-pos_targ(m))/180*2*pi)-1))+b_in2)*delta_t;
   end
   proba_in2 = proba_in2/((nu_targets+4)/6);



   
   
   for j=1:nu_trial   
      if round(j/50)==j/50
         fprintf('%d  ',j);
      end
       
      %%% pretrial       
           
      %% input spike trains during pretrial.
      %% 21 and 22: input from visual targets for training and testing set
      %% trial_dur is used here out of lazinness. It should be
      %% pre_trial_dur but we recycle the variable that will be used again
      %% during the trial.
      spikes_in21(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));
      spikes_in22(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));
 
      %%% proba_out01: probability of firing of output neuron during
      %%% pretrial in response to visual target alone. Set to proba_in2.
      %%% proba_out02: same for testing set
      proba_out01 = proba_in2*ones(1,pre_trial_dur)/delta_t;
      proba_out02 = proba_in2*ones(1,pre_trial_dur)/delta_t;
      
      %%% spikes_out01: pretrial output spike train for training set
      %%% spikes_out11: pretrial output spike train for testing set
      %%% during pretrial, output neurons simply respond like visual units
      spikes_out01 = spikes_in21(:,1:pre_trial_dur);
      spikes_out02 = spikes_in22(:,1:pre_trial_dur);      

      proba_out_av0 = proba_out_av0 + proba_out01;
      proba_out50_pre(count,:) = proba_out01(50,:);
      proba_out1_pre(count,:) = proba_out01(1,:);

%%%trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proba_out11: probability of firing nu_neurons x nu_time_steps. Training
% set
% spikes_out11: spike train of LIP neurons [nu_neurons x nu_time_steps] Training
% set. Determined from proba_out11
% proba_out12: probability of firing nu_neurons x nu_time_steps. Testing
% set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      aux1_proba1 = zeros(N_out,trial_dur+2);
      aux1_proba1(:,1)= proba_out01(:,pre_trial_dur);
      aux2_proba1 = zeros(N_out,trial_dur+2);
      
      aux1_proba2 = zeros(N_out,trial_dur+2);
      aux1_proba2(:,1)= proba_out02(:,pre_trial_dur);
      aux2_proba2 = zeros(N_out,trial_dur+2);

      proba_out11= zeros(N_out,trial_dur);
      proba_out12= zeros(N_out,trial_dur);
 
      proba_out11(:,1)= proba_out01(:,pre_trial_dur);
      proba_out12(:,1)= proba_out02(:,pre_trial_dur);
      spikes_out11(:,1)= spikes_out01(:,pre_trial_dur);
      spikes_out12(:,1)= spikes_out02(:,pre_trial_dur);     
      
      %% input spike trains
      %% 11 and 12: Motion input for training and testing set
      %% 21 and 22: input from visual targets for training and testing set
      spikes_in11(:,1:trial_dur)= rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));
      spikes_in12(:,1:trial_dur)= rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));

      spikes_in21(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));
      spikes_in22(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));


   
      proba_out_av1 = proba_out_av1 + proba_out11;
      proba_out50(count,:) = proba_out11(50,:);
      proba_out1(count,:) = proba_out11(1,:);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% training data. Spike counts and average probabilities.
      %% size nu_unit x 2 nu_trial (nu_trial for each of the 2 positions)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      spikes_in_trial11(:,count) = sum(spikes_in11')';

      proba_out_trial01(:,count) = sum(proba_out11');
      proba_out_trial11(:,count) = sum(proba_out11(:,1:wind_size)')/wind_size;
      proba_out_trial21(:,count) = sum(proba_out11(:,101:100+wind_size)')/wind_size;
      proba_out_trial31(:,count) = sum(proba_out11(:,201:200+wind_size)')/wind_size;
      proba_out_trial41(:,count) = sum(proba_out11(:,trial_dur-wind_size+1:trial_dur)')/wind_size;

      spikes_out_trial01(:,count) = sum(spikes_in11')';
      spikes_out_trial11(:,count) = sum(spikes_in11(:,1:wind_size)')';
      spikes_out_trial21(:,count) = sum(spikes_in11(:,1:100+wind_size)')';
      spikes_out_trial31(:,count) = sum(spikes_in11(:,1:200+wind_size)')';
      spikes_out_trial41(:,count) = sum(spikes_in11(:,1:trial_dur)')';
      
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% testing data. Spike counts and average probabilities.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      spikes_in_trial12(:,count) = sum(spikes_in12')';
 
      proba_out_trial02(:,count) = sum(proba_out12');
      proba_out_trial12(:,count) = sum(proba_out12(:,1:wind_size)')/wind_size;
      proba_out_trial22(:,count) = sum(proba_out12(:,101:100+wind_size)')/wind_size;
      proba_out_trial32(:,count) = sum(proba_out12(:,201:200+wind_size)')/wind_size;
      proba_out_trial42(:,count) = sum(proba_out12(:,trial_dur-wind_size+1:trial_dur)')/wind_size;
      
    
      spikes_out_trial02(:,count) = sum(spikes_in12')';
      spikes_out_trial12(:,count) = sum(spikes_in12(:,1:wind_size)')';
      spikes_out_trial22(:,count) = sum(spikes_in12(:,1:100+wind_size)')';
      spikes_out_trial32(:,count) = sum(spikes_in12(:,1:200+wind_size)')';
      spikes_out_trial42(:,count) = sum(spikes_in12(:,1:trial_dur)')';

      count = count+1;
   end
   fprintf('\n');   
end


%%%%%%%%%%%%%%%%%%%%%%%%
% stats
%%%%%%%%%%%%%%%%%%%%%%%%
mean_in  = mean(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
mean_in2 = mean(spikes_in_trial11(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);
var_in = var(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);

mean_out = mean(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);
mean_out2 = mean(spikes_out_trial01(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);

mean_proba  = mean(proba_out_trial01(:,1:nu_trial)');
mean_proba2 = mean(proba_out_trial01(:,1:nu_trial)');
mean_proba11 = mean(proba_out_trial11(:,1:nu_trial)');
mean_proba12 = mean(proba_out_trial12(:,1:nu_trial)');
mean_proba21 = mean(proba_out_trial21(:,1:nu_trial)');
mean_proba22 = mean(proba_out_trial22(:,1:nu_trial)');
mean_proba31 = mean(proba_out_trial31(:,1:nu_trial)');
mean_proba32 = mean(proba_out_trial32(:,1:nu_trial)');
mean_proba41 = mean(proba_out_trial41(:,1:nu_trial)');
mean_proba42 = mean(proba_out_trial42(:,1:nu_trial)');

var_out = var(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);
cov_in =  cov(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
cov_out =  cov(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);

difference = mean(spikes_out_trial01(:,1:nu_trial)')...
            -mean(spikes_out_trial01(:,nu_trial+1:nu_trial*2)');
       
proba_out_av0 = proba_out_av0/(2*nu_trial);            
proba_out_av1 = proba_out_av1/(2*nu_trial);            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info_flag==1
   % output info rate
   [info_rate_out info_rate_train_out w_ole_rate_out est_rate_angle perc_right_rate_test] = ...
        compute_info(proba_out_trial01, proba_out_trial02, step_test);

    % output info spikes
   [info_out info_train_out w_ole_out est_angle perc_right_test] = ...
        compute_info(spikes_out_trial01, spikes_out_trial02, step_test);

   %input info
   [info_in info_train_in w_ole_in est_angle_in perc_right_test_in] = ...
        compute_info(spikes_in_trial11, spikes_in_trial12, step_test);

   %normalize for duration of trial
   info_in = info_in/(trial_dur*delta_t);
   info_out= info_out/(trial_dur*delta_t);
   info_rate_out= info_rate_out/(trial_dur*delta_t);
   info_train_in = info_train_in/(trial_dur*delta_t);
   info_train_out= info_train_out/(trial_dur*delta_t);
   info_rate_train_out= info_rate_train_out/(trial_dur*delta_t);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info_flag==1

   [info_out1 info_train_out1 w_ole_out1 est_angle1 perc_right_test1] = ...
    compute_info(spikes_out_trial11, spikes_out_trial12, step_test);

   [info_out2 info_train_out2 w_ole_out2 est_angle2 perc_right_test2] = ...
    compute_info(spikes_out_trial21, spikes_out_trial22, step_test);

   [info_out3 info_train_out3 w_ole_out3 est_angle3 perc_right_test3] = ...
    compute_info(spikes_out_trial31, spikes_out_trial32, step_test);

   [info_out4 info_train_out4 w_ole_out4 est_angle4 perc_right_test4] = ...
    compute_info(spikes_out_trial41, spikes_out_trial42, step_test);
 

   %normalize for duration of windows
   info_out1 = info_out1/(wind_size*delta_t);
   info_out2 = info_out2/(wind_size*delta_t);
   info_out3 = info_out3/(wind_size*delta_t);
   info_out4 = info_out4/(wind_size*delta_t);
   info_train_out1 = info_train_out1/(wind_size*delta_t);
   info_train_out2 = info_train_out2/(wind_size*delta_t);
   info_train_out3 = info_train_out3/(wind_size*delta_t);
   info_train_out4 = info_train_out4/(wind_size*delta_t);
   
   
   if shuffle_flag ==1
      shuffle_data;
      
      [info_shuff_out1 info_train_shuff_out1 w_ole_shuff_out1 est_shuff_angle1 perc_right_shuff_test1] = ...
          compute_info(spikes_out_shuff_trial11, spikes_out_shuff_trial12, step_test);

      [info_shuff_out2 info_train_shuff_out2 w_ole_shuff_out2 est_shuff_angle2 perc_right_shuff_test2] = ...
          compute_info(spikes_out_shuff_trial21, spikes_out_shuff_trial22, step_test);

      [info_shuff_out3 info_train_shuff_out3 w_ole_shuff_out3 est_shuff_angle3 perc_right_shuff_test3] = ...
          compute_info(spikes_out_shuff_trial31, spikes_out_shuff_trial32, step_test);

      [info_shuff_out4 info_train_shuff_out4 w_ole_shuff_out4 est_shuff_angle4 perc_right_shuff_test4] = ...
          compute_info(spikes_out_shuff_trial41, spikes_out_shuff_trial42, step_test);
      
          %normalize for duration of windows
          info_shuff_out1 = info_shuff_out1/(wind_size*delta_t);
          info_shuff_out2 = info_shuff_out2/(wind_size*delta_t);
          info_shuff_out3 = info_shuff_out3/(wind_size*delta_t);
          info_shuff_out4 = info_shuff_out4/(wind_size*delta_t);
          info_train_shuff_out1 = info_train_shuff_out1/(wind_size*delta_t);
          info_train_shuff_out2 = info_train_shuff_out2/(wind_size*delta_t);
          info_train_shuff_out3 = info_train_shuff_out3/(wind_size*delta_t);
          info_train_shuff_out4 = info_train_shuff_out4/(wind_size*delta_t);

   end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate activation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out= mean_out;
in = mean_proba;

plot(mean_proba,2*log(1+exp(mean_proba/1000)),'-g');

parameters = fminsearch(@(param) dist2(param,in,out),[1000 2]);

alpha= parameters(1);
amp = parameters(2);

hold on
plot(mean_proba,mean_out,'o');
plot(mean_proba,amp*log(1+exp(mean_proba/alpha)),'-');
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute true info 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if info_flag ==1
   ac_in = max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1));
   ac_in_der = -max_rate_in * K_in * sin((pos_in'-pos)/180*2*pi).*...
            exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))/180*2*pi;

   info_in_true = sum(ac_in_der.^2./ac_in);

   true_cov_in =diag(ac_in);

   local_cov_out = diag(mean_out2);

% 1e-6 is added to ensure invertibility in the next line when aux_info is
% zero
   aux_info = w_oo*local_cov_out*w_oo' + diag(diag(ones(N_out)*1e-6));
   info_out_rate_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);

   D_mat = diag(amp/alpha*exp(mean_proba/alpha)./(1+exp(mean_proba/alpha)));
   D_mat = D_mat/max(max(D_mat));
   D_mat_inv = inv(D_mat);

aux_info = D_mat_inv*local_cov_out*D_mat_inv + diag(diag(ones(N_out)*1e-6));

%aux_info = local_cov_out + diag(diag(ones(N_out)*1e-6));
info_out_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate width of population activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=mean_in;
[parameters val] = fminsearch(@(param) dist(param, data),[max(mean_in); K_in; peak; b_in]);

amp_fit= parameters(1);
K_fit = parameters(2);
peak_fit = parameters(3);
offset = parameters(4)
fit_in = amp_fit * exp(K_fit*(cos((pos_in'-peak_fit)/180*2*pi)-1)) + offset;

hwhh_in = acos(1/K_fit*log((amp_fit-exp(-2*K_fit))/(2*amp_fit))+1)/(2*pi)*180;


data=mean_out;
fit_out_init = max(mean_out)*exp(K_in*(cos((pos_in'-peak)/180*2*pi)-1));
[parameters val] = fminsearch(@(param) dist(param, data),[max(mean_out); K_in; peak; b_in])

amp_fit= parameters(1);
K_fit  = parameters(2);
peak_fit = parameters(3);
offset = parameters(4)

fit_out = amp_fit*exp(K_fit*(cos((pos_in'-peak_fit)/180*2*pi)-1))+offset;

hwhh_out = acos(1/K_fit*log((amp_fit-exp(-2*K_fit))/(2*amp_fit))+1)/(2*pi)*180;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdf11 = exp(h_mat*spikes_out_trial11);
pdf11 = pdf11./(ones(N_out,1)*sum(pdf11));

pdf21 = exp(h_mat*spikes_out_trial21);
pdf21 = pdf21./(ones(N_out,1)*sum(pdf21));

pdf31 = exp(h_mat*spikes_out_trial31);
pdf31 = pdf31./(ones(N_out,1)*sum(pdf31));

pdf41 = exp(h_mat*spikes_out_trial41);
pdf41 = pdf41./(ones(N_out,1)*sum(pdf41));

pdf11_shuff = exp(h_mat*spikes_out_shuff_trial11);
pdf11_shuff = pdf11_shuff./(ones(N_out,1)*sum(pdf11_shuff));

pdf21_shuff = exp(h_mat*spikes_out_shuff_trial21);
pdf21_shuff = pdf21_shuff./(ones(N_out,1)*sum(pdf21_shuff));

pdf31_shuff = exp(h_mat*spikes_out_shuff_trial31);
pdf31_shuff = pdf31_shuff./(ones(N_out,1)*sum(pdf31_shuff));

pdf41_shuff = exp(h_mat*spikes_out_shuff_trial41);
pdf41_shuff = pdf41_shuff./(ones(N_out,1)*sum(pdf41_shuff));


%%% compute distribution over ML estimates
ML_info(1) = 0;

[aux_ML ML11] = max(pdf11);
pdf_ML11 = hist(ML11,[1:100]);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf_ML11),[max(pdf_ML11); K_in; peak; 0]);
K_fit = parameters(2);
ML_info(2) = K_fit;

[aux_ML ML21] = max(pdf21);
pdf_ML21 = hist(ML21,[1:100]);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf_ML21),[max(pdf_ML21); K_in; peak; 0]);
K_fit = parameters(2);
ML_info(3) = K_fit;

[aux_ML ML31] = max(pdf31);
pdf_ML31 = hist(ML31,[1:100]);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf_ML31),[max(pdf_ML31); K_in; peak; 0]);
K_fit = parameters(2);
ML_info(4) = K_fit;

[aux_ML ML41] = max(pdf41);
pdf_ML41 = hist(ML41,[1:100]);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf_ML41),[max(pdf_ML41); K_in; peak; 0]);
K_fit = parameters(2);
ML_info(5) = K_fit;


%%% compute average pdf over time
pdf11_aver = sum(pdf11')/sum(sum(pdf11'));
pdf21_aver = sum(pdf21')/sum(sum(pdf21'));
pdf31_aver = sum(pdf31')/sum(sum(pdf31'));
pdf41_aver = sum(pdf41')/sum(sum(pdf41'));

pdf11_shuff_aver = sum(pdf11_shuff')/sum(sum(pdf11_shuff'));
pdf21_shuff_aver = sum(pdf21_shuff')/sum(sum(pdf21_shuff'));
pdf31_shuff_aver = sum(pdf31_shuff')/sum(sum(pdf31_shuff'));
pdf41_shuff_aver = sum(pdf41_shuff')/sum(sum(pdf41_shuff'));


%%% compute log odds over time
log_odds(1)=0;
log_odds(2) = mean(log10(pdf11(50,:)./pdf11(100,:)));
log_odds(3) = mean(log10(pdf21(50,:)./pdf21(100,:)));
log_odds(4) = mean(log10(pdf31(50,:)./pdf31(100,:)));
log_odds(5) = mean(log10(pdf41(50,:)./pdf41(100,:)));

%%% Compute information over time by fitting von mises functions to the pdf
%%% and using the K parameter to estimate information
pdf_info(1) = 0;

pdf_mean(2) = sum(pdf11_aver.*pos_out);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf11_aver),[max(pdf11_aver); K_in; peak; 0]);
K_fit = parameters(2);
pdf_info(2) = K_fit;

pdf_mean(3) = sum(pdf21_aver.*pos_out);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf21_aver),[max(pdf21_aver); K_in; peak; 0]);
K_fit = parameters(2);
pdf_info(3) = K_fit;

pdf_mean(4) = sum(pdf31_aver.*pos_out);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf31_aver),[max(pdf31_aver); K_in; peak; 0]);
K_fit = parameters(2);
pdf_info(4) = K_fit;

pdf_mean(5) = sum(pdf41_aver.*pos_out);
[parameters val] =...
    fminsearch(@(param) dist(param, pdf41_aver),[max(pdf41_aver); K_in; peak; 0]);
K_fit = parameters(2);
pdf_info(5) = K_fit;


%%% compute pdf at decision time based on spike counts over the last 100ms 
pdf_decis = exp(h_mat*last_spikes);
pdf_decis = pdf_decis./(ones(N_out,1)*sum(pdf_decis));


%%% compute performance on each trial
[junk pos_max] = max(pdf_decis);
pos_max = pos_max/N_out*180;
bound_max = peak+step_test*(-0.5)+ 180/nu_targets/2;
bound_min = peak+step_test*(-0.5)- 180/nu_targets/2;
perf_trial = (bound_min<pos_max).*(pos_max<bound_max).*(RT'>0);


conf_level = log10(pdf_decis(50,:)./pdf_decis(100,:));

%%% compute overall performance
perf = sum(perf_trial)/(sum(RT>0));

RT_cor= RT.*(RT>0).*(perf_trial'>0);
RT_incor = RT.*(RT>0).*(perf_trial'<1);

mean_RT = sum(RT.*(RT>0).*(perf_trial'>0))/(sum((RT>0).*(perf_trial'>0)));

mean_RT_incor = sum(RT.*(RT>0).*(perf_trial'<1))/(sum((RT>0).*(perf_trial'<1)));



%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%
if info_flag==1
fprintf('\n');
fprintf('True Info In test       : %2.3f\n',info_in_true);
fprintf('True Info rate Out test : %2.3f\n',info_out_rate_true);
fprintf('True Info Out test      : %2.3f\n',info_out_true);
fprintf('\n');

fprintf('Info In test           : %2.3f\n',info_in);
fprintf('Info In train          : %2.3f\n',info_train_in);
fprintf('Percentage classification correct: %2.2f\n',perc_right_test_in);
fprintf('\n');

fprintf('Info Out test          : %2.3f\n',info_out);
fprintf('Info Out train         : %2.3f\n',info_train_out);
fprintf('Percentage classification correct: %2.2f\n\n',perc_right_test);

fprintf('Info Out test 1   1-%d ms   : %2.3f  %2.3f\n',wind_size,info_out1,info_train_out1);
fprintf('Info Out test 2 101-%d ms   : %2.3f  %2.3f\n',100+wind_size,info_out2,info_train_out2);
fprintf('Info Out test 3 201-%d ms  : %2.3f  %2.3f\n',200+wind_size,info_out3,info_train_out3);
fprintf('Info Out test 4 %d-500 ms  : %2.3f  %2.3f\n',501-wind_size,info_out4,info_train_out4);
fprintf('\n');


fprintf('Info rate Out test          : %2.3f\n',info_rate_out);
fprintf('Info rate Out train         : %2.3f\n',info_rate_train_out);
fprintf('Percentage classification correct: %2.2f\n',perc_right_rate_test);
fprintf('\n');

end 

fprintf('Max firing rate: %2.3f\n',max(mean_out));
fprintf('\n');

fprintf('Norm of W_oo : %2.6f\n',w_oo(:,1)'*w_oo(:,1));
fprintf('Norm of W_oi: %2.6f\n',w_oi(:,1)'*w_oi(:,1));
fprintf('\n');
fprintf('\n');

fprintf('Width in      : %2.3f\n',hwhh_in);
fprintf('Width out      : %2.3f\n',hwhh_out);
fprintf('\n');

fprintf('Performance: %f\n', perf);
fprintf('Mean RT correct  : %f ms\n', mean_RT);
fprintf('Mean RT incorrect: %f ms\n', mean_RT_incor);




subplot(221)


if nu_targets==2
   lik1 = [1 pdf11_aver(50) pdf21_aver(50) pdf31_aver(50) pdf41_aver(50)];
   lik2 = [1 pdf11_aver(100) pdf21_aver(100) pdf31_aver(100) pdf41_aver(100)];
   loglik= log10(lik1./lik2);
   lik1_shuff = [1 pdf11_shuff_aver(50) pdf21_shuff_aver(50) pdf31_shuff_aver(50) pdf41_shuff_aver(50)];
   lik2_shuff = [1 pdf11_shuff_aver(100) pdf21_shuff_aver(100) pdf31_shuff_aver(100) pdf41_shuff_aver(100)];
   loglik_shuff= log10(lik1_shuff./lik2_shuff);

   plot([1 100 200 300 trial_dur], loglik,'*-');
   hold on
   plot([1 100 200 300 trial_dur], loglik_shuff,'r*-'); 
   hold off
   title('Log odds');
elseif info_flag==1
   plot([1 100 200 300 trial_dur-1], [0 info_out1 info_out2 info_out3 info_out4],'-o')
   hold on
   plot([1 100 200 300 trial_dur-1], [0 info_shuff_out1 info_shuff_out2 info_shuff_out3 info_shuff_out4],'r-o')
   hold off
   title('Fisher information');
end




subplot(222)
% plot(mean_in,var_in,'.');
% hold on
% plot(mean_out,var_out,'.r');
% plot([1 max(mean_out)],[1 max(mean_out)])
% hold off

plot(pdf11(:,1),'.r');
hold on 
plot(pdf21(:,1),'.b');
plot(pdf31(:,1),'.g');
%plot(pdf41(:,1),'.c');
plot(pdf11_aver,'r')
plot(pdf21_aver,'b')
plot(pdf31_aver,'g')
%plot(pdf41_aver,'c')

%errorbar([1:N_out],sum(pdf31')/sum(sum(pdf31')),std(pdf31'),'r')
%errorbar([1:N_out],sum(pdf41')/sum(sum(pdf41')),std(pdf41'),'c')
hold off
axis tight
ylabel('P(s|r)')
xlabel('Stimulus')



% plot(mean_proba11(:,1),'.r');
% hold on 
% plot(pdf21(:,1),'.b');
% plot(pdf31(:,1),'.g');
% plot(pdf41(:,1),'.c');
% plot(sum(pdf11')/sum(sum(pdf11')),'r')
% plot(sum(pdf21')/sum(sum(pdf21')),'b')
% plot(sum(pdf31')/sum(sum(pdf31')),'g')
% plot(sum(pdf41')/sum(sum(pdf41')),'c')
% hold off

% aux1111 = log(sum(pdf11')/sum(sum(pdf11')))
% aux1111 = aux1111-min(aux1111);
% aux1111 = aux1111/max(aux1111);
% aux1112 = log(sum(pdf21')/sum(sum(pdf21')))
% aux1112 = aux1112-min(aux1112);
% aux1112 = aux1112/max(aux1112);
% aux1113 = log(sum(pdf31')/sum(sum(pdf31')))
% aux1113 = aux1113-min(aux1113);
% aux1113 = aux1113/max(aux1113);
% aux1114 = log(sum(pdf41')/sum(sum(pdf41')))
% aux1114 = aux1114-min(aux1114);
% aux1114 = aux1114/max(aux1114);
% 
% plot(aux1111,'r')
% hold on 
% plot(aux1112,'b')
% plot(aux1113,'g')
% plot(aux1114,'c')
% hold off





subplot(223)
plot(pos_in,mean_in,'o')
hold on
plot(pos_out,mean_out,'or')
plot(pos_in,fit_in,'-b')
plot(pos_in,fit_out,'-r')
hold off
max_plot=max(max(mean_out),max(mean_in));
axis([0 180 0 max_plot*1.2]);

plot(mean(spikes_out_trial01),'or')


aux111 = mean(spikes_out_trial11');
aux111 = aux111-min(aux111);
aux111 = aux111/max(aux111);
aux112 = mean(spikes_out_trial21');
aux112 = aux112-min(aux112);
aux112 = aux112/max(aux112);
aux113 = mean(spikes_out_trial31');
aux113 = aux113-min(aux113);
aux113 = aux113/max(aux113);
aux114 = mean(spikes_out_trial41');
aux114 = aux114-min(aux114);
aux114 = aux114/max(aux114);

plot(aux111,'or')
hold on
plot(aux112,'ob')
plot(aux113,'og')
plot(aux114,'oc')
hold off


% plot([0 100 200 300 400],log_odds,'o-');
% xlabel('Time (ms)')
% ylabel('Log odds')




subplot(224)
%surf([0:180/(N_out-1):180], [0:180/(N_out-1):180],cov_out-diag(diag(cov_out)))
%shading interp
aux_proba= [proba_out_av0(:,pre_trial_dur-49:10:pre_trial_dur) proba_out_av1(:,1:10:trial_dur)];
%surfl(aux_proba)

plot([-49:10:trial_dur],aux_proba(50,:))
hold on
plot([-49:10:trial_dur],aux_proba(100,:),'g')
plot([-49:10:trial_dur],aux_proba(25,:),'c')
hold off
axis([-49 trial_dur 0 60] )
title('Conditioned on neuron')
% 
% proba_out50_cor = proba_out50.*(perf_trial'*ones(1,trial_dur));
% proba_out1_cor = proba_out1.*(perf_trial'*ones(1,trial_dur));
% 
% neuron_traj = sum(proba_out50_cor)/(sum(perf_trial));
% antineuron_traj = sum(proba_out1_cor)/(sum(perf_trial));
% 
% plot(neuron_traj)
% hold on
% plot(antineuron_traj,'g');
% hold off
% %axis([0 400 20 40])
% ylabel('Activity (Hz)');
% xlabel('Time (ms)');
% title('Conditioned on correct')

save datalast RT conf_level perf_trial perf mean_RT mean_RT_incor coherence neuron_traj antineuron_traj

beep
toc

