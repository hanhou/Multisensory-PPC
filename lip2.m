%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
%% this version is used for batch simulations
%% to obtain the RT and performance curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load h

parameters;


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
   pos = peak+step_test*(l-step_test/2);
   fprintf('%d  pos: %4.2f\n  ',l, pos);
   
%    %%% proba_in: proba of firing of input neurons in response to motion
     max_rate_in = r_spont_MT + b_pref*coherence;
     b_in = r_spont_MT + b_null*coherence;
     proba_in = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
     
     max_rate_in = r_spont_MT + b_pref*coherence*2;
     b_in = r_spont_MT + b_null*coherence*2;
     proba_in_temp = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 

%      max_rate_in = b_pref*coherence;
%      proba_in = (max_rate_in/2*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t; 
   
%     proba_in_temp = (max_rate_in/2*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t; ; 
%     proba_in_temp2 = (max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t; ; 
 
%    max_rate_in = r_spont_MT + b_pref*coherence;
%    b_in = r_spont_MT + b_null*coherence;
%    proba_in_temp = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
   
   %%% proba_in2: proba of firing of input neurons in response to visual
   %%% targets
   pos_targ = pos+[0:180/nu_targets:179];
   proba_in2 = zeros(N_out,1);
   for m=1:nu_targets
       aux(:,m)= ((max_rate_in2-b_in2)*...
           exp(K_in2*(cos((pos_in'-pos_targ(m))/180*2*pi)-1))+b_in2)*delta_t;
       proba_in2 = proba_in2 + ((max_rate_in2-b_in2)*...
           exp(K_in2*(cos((pos_in'-pos_targ(m))/180*2*pi)-1))+b_in2)*delta_t;
   end
   proba_in2 = proba_in2/(slope_norm*(nu_targets/2)+dc_norm);

   
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
% proba_out{1}: probability of firing nu_neurons x nu_time_steps. Training
% set
% spikes_out{1}: spike train of LIP neurons [nu_neurons x nu_time_steps] Training
% set. Computed from proba_out{1}.
% proba_out{2}: probability of firing nu_neurons x nu_time_steps. Testing
% set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      aux1_proba1 = zeros(N_out,trial_dur+2);
      aux1_proba1(:,1)= proba_out01(:,pre_trial_dur);
      aux2_proba1 = zeros(N_out,trial_dur+2);
      
      aux1_proba2 = zeros(N_out,trial_dur+2);
      aux1_proba2(:,1)= proba_out02(:,pre_trial_dur);
      aux2_proba2 = zeros(N_out,trial_dur+2);

      proba_out{1} = zeros(N_out,trial_dur);
      proba_out{2} = zeros(N_out,trial_dur);
 
      proba_out{1}(:,1)= proba_out01(:,pre_trial_dur);
      proba_out{2}(:,1)= proba_out02(:,pre_trial_dur);

      spikes_out{1}(:,:)= zeros(N_out,trial_dur);
      spikes_out{2}(:,:)= zeros(N_out,trial_dur);
      spikes_out{1}(:,1)= spikes_out01(:,pre_trial_dur);
      spikes_out{2}(:,1)= spikes_out02(:,pre_trial_dur);     
      
      %% spike_in: input spike trains due to motion. Two sets for 
      %% training and testing sets.
      for m=1:2
         aux_proba_in = proba_in + w_mt*randn(N_in,1);
         spikes_in{m}(:,1:trial_dur)= rand(N_in,trial_dur)<(aux_proba_in*ones(1,trial_dur));    
%          aux_proba_in = proba_in_temp + w_mt*randn(N_in,1);
%          spikes_in{m}(:,101:200)= rand(N_in,100)<(aux_proba_in*ones(1,100));
      end

      spikes_in21(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));
      spikes_in22(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));

      count_time50=1;
      count_time100=1;
      k=1;
      while k<=trial_dur-1        
         aux1_proba1(:,k+1) = (1-delta_t/time_const_in)*aux1_proba1(:,k)...
                              +1/time_const_in*(w_oi*spikes_in{1}(:,k)...
                                        +w_oi2*att_gain2*spikes_in21(:,k));
         aux2_proba1(:,k+1) = (1-delta_t/time_const_out)*aux2_proba1(:,k)+...
                              +1/time_const_out*(w_oo*spikes_out{1}(:,k));                         
         proba_out{1}(:,k+1) = aux1_proba1(:,k+1) + aux2_proba1(:,k+1) + b_out;
                            
         aux1_proba2(:,k+1) = (1-delta_t/time_const_in)*aux1_proba2(:,k)...
                              +1/time_const_in*(w_oi*spikes_in{2}(:,k)...
                                        +w_oi2*att_gain2*spikes_in22(:,k));
         aux2_proba2(:,k+1) = (1-delta_t/time_const_out)*aux2_proba2(:,k)...
                              +1/time_const_out*(w_oo*spikes_out{2}(:,k));
         proba_out{2}(:,k+1) = aux1_proba2(:,k+1) + aux2_proba2(:,k+1) + b_out;
                    
         spikes_out{1}(:,k+1)  = (rand(N_out,1) < delta_t*(proba_out{1}(:,k+1)-threshold));
         spikes_out{2}(:,k+1) = (rand(N_out,1) < delta_t*(proba_out{2}(:,k+1)-threshold));

         %% variable used for stopping the integration
%          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
%                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));    
         decision_ac(:,k+1) = proba_out{1}(:,k+1);
                              
         % collect spike times into a cell
         if spikes_out{1}(50,k+1) == 1
             spike_time50{count}(count_time50) = k+1;
             count_time50 = count_time50+1;
         end
               
         if spikes_out{1}(100,k+1) == 1
             spike_time100{count}(count_time100) = k+1;
             count_time100 = count_time100+1;
         end
        
         %termination 
         if ((max(decis_flag*decision_ac(:,k+1))> decis_thres) && (RT(count)==0)) || (k==trial_dur-1)
             RT(count,1) = k;
             if k>wind_size-1
                last_spikes(:,count) = sum(spikes_out{1}(:,k-wind_size+1:k)')';
             else
                last_spikes(:,count) = sum(spikes_out{1}(:,1:k)')';
             end
             last_proba(:,count) = proba_out{1}(:,k);
             k=trial_dur;
         end
         k=k+1;
      end 
      
       
      proba_out_av1 = proba_out_av1 + proba_out{1};
      proba_out50(count,:) = proba_out{1}(50,:);
      proba_out75(count,:) = proba_out{1}(75,:);
      proba_out100(count,:) = proba_out{1}(100,:);
      spikes_out50(count,:) = spikes_out{1}(50,:);
      spikes_out75(count,:) = spikes_out{1}(75,:);
      spikes_out100(count,:) = spikes_out{1}(100,:);
      spikes_in50(count,:) = spikes_in{1}(50,:);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% training and testing data. Spike counts and average probabilities.
      %% size nu_unit x 2 nu_trial (nu_trial for each of the 2 positions)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      spikes_in_trial11(:,count) = sum(spikes_in{1}')';
      spikes_in_trial12(:,count) = sum(spikes_in{2}')';
       
      spikes_out_trial01(:,count) = sum(spikes_out{1}')';
      spikes_out_trial02(:,count) = sum(spikes_out{2}')';
      
      proba_out_trial01(:,count) = sum(proba_out{1}(:,1:trial_dur)')'/wind_size;         
      proba_out_trial02(:,count) = sum(proba_out{2}(:,1:trial_dur)')'/wind_size;         
 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% collect count in small time windows
      %% m=1 training data, m=2 testing data
      %% n?: index of time window
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      for m=1:2
      for n=1:nu_window
         proba_count{m,n}(:,count)= ...
             sum(proba_out{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')'/wind_size;         
         spikes_count{m,n}(:,count)= ...
             sum(spikes_out{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')';
         spikes_count_in{m,n}(:,count)= ...
             sum(spikes_in{m}(:,1:wind_offset+wind_size*n)')';
      end
      end
 
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

for m=1:2
for n=1:nu_window
   mean_proba_count{m,n}=  mean(proba_count{m,n}(:,1:nu_trial)');
end
end

var_out = var(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);
cov_in =  cov(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
corr_in =  cov_in./(diag(cov_in).^0.5*(diag(cov_in)').^0.5);

%compute the average correlation coefficients in MT
aux_mask = ones(N_in,1)*[1:N_in];
mask = ((1-cos(abs(aux_mask-aux_mask')/N_in*2*pi))/2)<.5;
mask = mask.*(1-eye(N_in));
cor90 = sum(sum(corr_in.*mask))/sum(sum(mask));
mask = ((1-cos(abs(aux_mask-aux_mask')/N_in*2*pi))/2)>.5;
cor180 = sum(sum(corr_in.*mask))/sum(sum(mask));

cov_out =  cov(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);

difference = mean(spikes_out_trial01(:,1:nu_trial)')...
            -mean(spikes_out_trial01(:,nu_trial+1:nu_trial*2)');
       
proba_out_av0 = proba_out_av0/(2*nu_trial);            
proba_out_av1 = proba_out_av1/(2*nu_trial);            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info_in_flag==1
   %input info
   [info_in info_train_in w_ole_in est_angle_in perc_right_test_in] = ...
        compute_info(spikes_in_trial11, spikes_in_trial12, step_test);

   %normalize for duration of trial
   info_in = info_in/(trial_dur*delta_t);
   info_train_in = info_train_in/(trial_dur*delta_t);
end

if info_out_flag==1
    % output info spikes
   [info_out info_train_out w_ole_out est_angle perc_right_test] = ...
        compute_info(spikes_out_trial01, spikes_out_trial02, step_test);

   %normalize for duration of trial
   info_out= info_out/(trial_dur*delta_t);
   info_train_out= info_train_out/(trial_dur*delta_t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_outx=zeros(1,nu_window);
w_ole_outx = zeros(N_out+1,nu_window);
est_anglex= zeros(nu_window,nu_trial*2);

%%% information in MT
if info_in_flag==1
   for m=1:nu_window
      [info_inx(m) info_train_inx(m) w_ole_inx(:,m) est_angle_inx(m,:) perc_right_test_inx(m)] = ...
         compute_info(spikes_count_in{1,m}, spikes_count_in{2,m}, step_test);
       %normalize for duration of window
      info_inx(m) = info_inx(m) /(wind_size*delta_t);
      info_train_inx(m) = info_train_inx(m)/(wind_size*delta_t);      
   end
end   
  
%%%information in LIP
if info_out_flag==1
   for m=1:nu_window
      [info_outx(m) info_train_outx(m) w_ole_outx(:,m) est_anglex(m,:) perc_right_testx(m)] = ...
         compute_info(spikes_count{1,m}, spikes_count{2,m}, step_test);
       %normalize for duration of window
      info_outx(m) = info_outx(m) /(wind_size*delta_t);
      info_train_outx(m) = info_train_outx(m)/(wind_size*delta_t);
      
      [info_proba(m) info_train_proba(m) w_ole_proba(:,m) est_angle_proba(m,:) perc_right_test_proba(m)] = ...
         compute_info(proba_count{1,m}, proba_count{2,m}, step_test);
       %normalize for duration of window
      info_proba(m) = info_proba(m) /(wind_size*delta_t);
      info_train_proba(m) = info_train_proba(m)/(wind_size*delta_t);
   end
   
   
   
   if shuffle_flag ==1
      for j=1:N_out
         for m=1:2
         for n=1:nu_window
            spikes_count_shuff{m,n}(j,1:nu_trial) = ...
                spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
            spikes_count_shuff{m,n}(j,nu_trial+1:nu_trial*2) = ...
                spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         end
         end
      end
      
      for m=1:nu_window
         [info_shuff_out(m) info_train_shuff_out(m)...
             w_ole_shuff_out est_shuff_angle perc_right_shuff_test(m)] = ...
             compute_info(spikes_count_shuff{1,m}, spikes_count_shuff{2,m}, step_test);
         %normalize for duration of window
         info_shuff_out(m) = info_shuff_out(m)/(wind_size*delta_t);
         info_train_shuff_out(m) = info_train_shuff_out(m)/(wind_size*delta_t);
      end
   end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows using beck's grad descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if info_out_flag==1
%    for m=1:nu_window
%       [binfo_train_out(m) binfo_out(m) bw_ole_out(:,m) bperc_right_train(m) bperc_right_test(m)] = ...
%          infobeck(spikes_count{1,m}, spikes_count{2,m}, step_test, 0.01);
%        %normalize for duration of window
%       binfo_out(m) = binfo_out(m) /(wind_size*delta_t);
%       binfo_train_out(m) = binfo_train_out(m)/(wind_size*delta_t);
%       
%       [binfo_train_proba(m) binfo_proba(m) bw_ole_proba(:,m) bperc_right_train_proba(m) bperc_right_test_proba(m)] = ...
%          infobeck(proba_count{1,m}, proba_count{2,m}, step_test, 0.01);
%        %normalize for duration of window
%       binfo_proba(m) = binfo_proba(m) /(wind_size*delta_t);
%       binfo_train_proba(m) = binfo_train_proba(m)/(wind_size*delta_t);
%    end
% end

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

if info_out_flag ==1
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
fprintf('Estimate the width of the population activity \n');
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
[parameters val] = fminsearch(@(param) dist(param, data),[max(mean_out); K_in; peak; b_in]);

amp_fit= parameters(1);
K_fit  = parameters(2);
peak_fit = parameters(3);
offset = parameters(4)

fit_out = amp_fit*exp(K_fit*(cos((pos_in'-peak_fit)/180*2*pi)-1))+offset;

hwhh_out = acos(1/K_fit*log((amp_fit-exp(-2*K_fit))/(2*amp_fit))+1)/(2*pi)*180;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Estimate the pdfs \n');

for m=1:nu_window
   pdf{m} = exp(h_mat*spikes_count{1,m});
   pdf{m} = pdf{m}./(ones(N_out,1)*sum(pdf{m}));
   pdf_aver(m,:) = sum(pdf{m}')/sum(sum(pdf{m}'));
end

if shuffle_flag==1
    for m=1:nu_window 
       pdf_shuff{m} = exp(h_mat*spikes_count_shuff{1,m});
       pdf_shuff{m} = pdf_shuff{m}./(ones(N_out,1)*sum(pdf_shuff{m}));
       pdf_shuff_aver{m} = sum(pdf_shuff{m}')/sum(sum(pdf_shuff{m}'));
    end
end

%%% compute distribution over ML estimates
ML_info(1) = 0;

for m=1:nu_window
   [aux_ML ML{m}] = max(pdf{m});
   pdf_ML{m} = hist(ML{m},[1:100]);
   [parameters val] =...
    fminsearch(@(param) dist(param, pdf_ML{m}),[max(pdf_ML{m}); K_in; peak; ,min(pdf_ML{m})]);
   K_fit = parameters(2);
   vardeg = 1/K_fit*(180/(2*pi))^2;
   ML_info(m+1) = 1/vardeg;
   ML_info(m+1) =ML_info(m+1)/(wind_size*delta_t);
end

%%% compute log odds over time
log_odds(1)=0;
for m=1:nu_window
    log_odds(m+1) = mean(log10(pdf{m}(50,:)./pdf{m}(100,:)));
end


%%% Compute information over time by fitting von mises functions to the pdf
%%% and using the K parameter to estimate information
pdf_info(1) = 0;
K_fit(1) = 10^-30;

for m=1:nu_window
   pdf_mean(m+1) = sum(pdf_aver(m,:).*pos_out);
   
   [parameters val] =...
    fminsearch(@(param) dist(param, pdf_aver(m,:)),[max(pdf_aver(m,:));...
               K_in; peak; min((pdf_aver(m,:)))]);
   K_fit(m+1) = parameters(2);
   vardeg = 1/K_fit(m+1)*(180/(2*pi))^2;
   pdf_est(m,:)= parameters(1)*exp(K_fit(m+1)*...
         (cos((pos_out-parameters(3))/180*2*pi)-1))+parameters(4);
   pdf_info(m+1) = 1/vardeg;
   pdf_info(m+1) = pdf_info(m+1)/(wind_size*delta_t);
end



%%% compute pdf at decision time based on spike counts over the last 100ms 
pdf_decis = exp(h_mat*last_spikes);
pdf_decis = pdf_decis./(ones(N_out,1)*sum(pdf_decis));


%%% compute performance on each trial
%[junk pos_max] = max(pdf_decis);
[junk pos_max] = max(last_proba);
pos_max = pos_max/N_out*180;
bound_max = peak+step_test*(-0.5)+ 180/nu_targets/2;
bound_min = peak+step_test*(-0.5)- 180/nu_targets/2;
perf_trial = (bound_min<pos_max).*(pos_max<bound_max).*(RT'>0);


conf_level = log10(pdf_decis(50,:)./pdf_decis(100,:));
conf_level_proba = last_proba(50,:)-last_proba(100,:);

%%% compute overall performance
perf = sum(perf_trial)/(sum(RT>0));

RT_cor= RT.*(RT>0).*(perf_trial'>0);
RT_incor = RT.*(RT>0).*(perf_trial'<1);

mean_RT = sum(RT.*(RT>0).*(perf_trial'>0))/(sum((RT>0).*(perf_trial'>0)));

mean_RT_incor = sum(RT.*(RT>0).*(perf_trial'<1))/(sum((RT>0).*(perf_trial'<1)));



%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%
if info_in_flag==1
    fprintf('\n');
    fprintf('Info In test           : %2.3f\n',info_in);
    fprintf('Info In train          : %2.3f\n',info_train_in);
end 

if info_out_flag==1
    fprintf('True Info rate Out test : %2.3f\n',info_out_rate_true);
    fprintf('True Info Out test      : %2.3f\n',info_out_true);
    fprintf('\n');

    fprintf('Percentage classification correct: %2.2f\n',perc_right_test_in);
    fprintf('\n');

    fprintf('Info Out test          : %2.3f\n',info_out);
    fprintf('Info Out train         : %2.3f\n',info_train_out);
    fprintf('Percentage classification correct: %2.2f\n\n',perc_right_test);
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


plot_results;