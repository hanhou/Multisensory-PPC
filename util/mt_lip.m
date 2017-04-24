%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
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
step_out = 180/N_out;
pos_out=[0:step_out:179.9999];
nu_targets = 2;
decis_thres = 28;

peak = 90;
pre_trial_dur= 50; %% must be 50 or more
trial_dur = 401;

if pre_trial_dur>trial_dur
   error('pre_trial_dur cannot be longer than trial_dur'); 
end

delta_t = 1e-3;
% difference in s between the two test stimuli
step_test =6.4;

nu_trial = 1000;
coherence = 20;

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

%%% duration of the window used to compute Fisher
wind_size = 100;
nu_window = 4;




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
   max_rate_in = r_spont_MT + b_pref*coherence;
   b_in = r_spont_MT + b_null*coherence;
   proba_in = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
   
%    max_rate_in = r_spont_MT + b_pref*coherence*2;
%    b_in = r_spont_MT + b_null*coherence*2;
%    proba_in_aux1 = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
% 
%    max_rate_in = r_spont_MT + b_pref*coherence*4;
%    b_in = r_spont_MT + b_null*coherence*4;
%    proba_in_aux2 = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
%    
%    max_rate_in = r_spont_MT + b_pref*coherence*8;
%    b_in = r_spont_MT + b_null*coherence*8;
%    proba_in_aux3 = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t; 
   
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
   proba_in2 = proba_in2/((nu_targets)/2);


   
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
% set. Determined from proba_out{1}
% proba_out{2}: probability of firing nu_neurons x nu_time_steps. Testing
% set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      aux1_proba1 = zeros(N_out,trial_dur+2);
      aux1_proba1(:,1)= proba_out01(:,pre_trial_dur);
      aux2_proba1 = zeros(N_out,trial_dur+2);
      
      aux1_proba2 = zeros(N_out,trial_dur+2);
      aux1_proba2(:,1)= proba_out02(:,pre_trial_dur);
      aux2_proba2 = zeros(N_out,trial_dur+2);

      proba_out{1}= zeros(N_out,trial_dur);
      proba_out{2}= zeros(N_out,trial_dur);
 
      proba_out{1}(:,1)= proba_out01(:,pre_trial_dur);
      proba_out{2}(:,1)= proba_out02(:,pre_trial_dur);
      spikes_out{1}(:,1)= spikes_out01(:,pre_trial_dur);
      spikes_out{2}(:,1)= spikes_out02(:,pre_trial_dur);     
      
      %% input spike trains
      %% 11 and 12: Motion input for training and testing set
      %% 21 and 22: input from visual targets for training and testing set
      spikes_in11(:,1:trial_dur)=   rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));
      spikes_in12(:,1:trial_dur)=   rand(N_in,trial_dur)<(proba_in*ones(1,trial_dur));

      spikes_in21(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));
      spikes_in22(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2*ones(1,trial_dur));

      count_time50=1;
      count_time100=1;

      
       
      proba_out_av1 = proba_out_av1 + proba_out{1};
      proba_out50(count,:) = proba_out{1}(50,:);
      proba_out1(count,:) = proba_out{1}(1,:);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% training and testing data. Spike counts and average probabilities.
      %% size nu_unit x 2 nu_trial (nu_trial for each of the 2 positions)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      spikes_in_trial11(:,count) = sum(spikes_in11')';
      spikes_in_trial12(:,count) = sum(spikes_in12')';
       
      spikes_out_trial01(:,count) = sum(spikes_out{1}')';
      spikes_out_trial02(:,count) = sum(spikes_out{2}')';
      
      proba_out_trial01(:,count) = sum(proba_out{1}(:,1:trial_dur)')'/wind_size;         
      proba_out_trial02(:,count) = sum(proba_out{2}(:,1:trial_dur)')'/wind_size;         
 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% collect count in small time windows
      %% m=1 training data, m=2 testing data
      %% n?: index of time window
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      for n=1:nu_window
         proba_count{1,n}(:,count)= ...
             sum(spikes_in11(:,1:wind_size*n)')';
         spikes_count{1,n}(:,count)= ...
             sum(spikes_in11(:,1:wind_size*n)')';
         proba_count{2,n}(:,count)= ...
             sum(spikes_in12(:,1:wind_size*n)')';
         spikes_count{2,n}(:,count)= ...
             sum(spikes_in12(:,1:wind_size*n)')';
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
cov_out =  cov(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);

difference = mean(spikes_out_trial01(:,1:nu_trial)')...
            -mean(spikes_out_trial01(:,nu_trial+1:nu_trial*2)');
       
proba_out_av0 = proba_out_av0/(2*nu_trial);            
proba_out_av1 = proba_out_av1/(2*nu_trial);            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info_flag==1
    % output info spikes
   [info_out info_train_out w_ole_out est_angle perc_right_test] = ...
        compute_info(spikes_out_trial01, spikes_out_trial02, step_test);

   %input info
   [info_in info_train_in w_ole_in est_angle_in perc_right_test_in] = ...
        compute_info(spikes_in_trial11, spikes_in_trial12, step_test);

   %normalize for duration of trial
   info_in = info_in/(trial_dur*delta_t);
   info_out= info_out/(trial_dur*delta_t);
   info_train_in = info_train_in/(trial_dur*delta_t);
   info_train_out= info_train_out/(trial_dur*delta_t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_outx=zeros(1,nu_window);
w_ole_outx = zeros(N_out+1,nu_window);
est_anglex= zeros(nu_window,nu_trial*2);

if info_flag==1
   for m=1:nu_window
      [info_outx(m) info_train_outx(m) w_ole_outx(:,m) est_anglex(m,:) perc_right_testx(m)] = ...
         compute_info(spikes_count{1,m}, spikes_count{2,m}, step_test);
       %normalize for duration of window
      info_outx(m) = info_outx(m) /(wind_size*delta_t);
      info_train_outx(m) = info_train_outx(m)/(wind_size*delta_t);
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

for m=1:4
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
   lik1 = [1 pdf_aver(1,50) pdf_aver(2,50) pdf_aver(3,50) pdf_aver(4,50)];
   lik2 = [1 pdf_aver(1,100) pdf_aver(2,100) pdf_aver(3,100) pdf_aver(4,100)];
   loglik= log10(lik1./lik2);
%    if shuffle_flag ==1
%       lik1_shuff = [1 pdf_shuff_aver{1}(50) pdf_shuff_aver{2}(50) pdf_shuff_aver{3}(50) pdf_shuff_aver{4}(50)];
%       lik2_shuff = [1 pdf11_shuff_aver{1}(100) pdf_shuff_aver{2}(100) pdf_shuff_aver{4}(100) pdf_shuff_aver{4}(100)];
%       loglik_shuff= log10(lik1_shuff./lik2_shuff);
%    end
   
   plot([1 100 200 300 trial_dur], loglik,'*-');
%    if shuffle_flag==1
%        hold on
%        plot([1 100 200 300 trial_dur], loglik_shuff,'r*-');
%        hold off
%    end
   title('Log odds');
elseif info_flag==1 
   plot([1 100 200 300 400],[0 info_outx],'r-o')
   hold on
   if shuffle_flag==1
      plot([1 100 200 300 trial_dur-1], [0 info_shuff_out],'b-*')
   end
   plot([1 100 200 300 400],[0 info_train_outx],'g:o')
   hold off
   title('Fisher information');
end




subplot(222)
plot(pdf{1}(:,1),'.r');
hold on 
plot(pdf{2}(:,1),'.b');
plot(pdf{3}(:,1),'.g');
plot(pdf{4}(:,1),'.c');
plot(pdf_aver(1,:),'r')
plot(pdf_aver(2,:),'b')
plot(pdf_aver(3,:),'g')
plot(pdf_aver(4,:),'c')
hold off

axis tight
ylabel('P(s|r)')
xlabel('Stimulus')





subplot(223)
plot(pos_in,mean_in,'o')
hold on
plot(pos_out,mean_out,'or')
plot(pos_in,fit_in,'-b')
plot(pos_in,fit_out,'-r')
hold off
max_plot=max(max(mean_out),max(mean_in));
axis([0 180 0 max_plot*1.2]);


for m=1:nu_window
    aux11{m} = mean(proba_count{1,m}');
%     aux11{m} = aux11{m}-min(aux11{m});
%     aux11{m} = aux11{m}/max(aux11{m});
end

plot(aux11{1},'or')
hold on
plot(aux11{2},'ob')
plot(aux11{3},'og')
plot(aux11{4},'oc')
hold off




subplot(224)
%surf([0:180/(N_out-1):180], [0:180/(N_out-1):180],cov_out-diag(diag(cov_out)))
%shading interp
aux_proba= [proba_out_av0(:,pre_trial_dur-49:10:pre_trial_dur) proba_out_av1(:,1:10:trial_dur)];
%surfl(aux_proba)

plot([-49:10:trial_dur],aux_proba(50,:))
hold on
plot([-49:10:trial_dur],aux_proba(100,:),'g')
plot([-49:10:trial_dur],aux_proba(25,:),'c')
plot([-49:10:trial_dur],aux_proba(50,:)-(aux_proba(100,:).*(aux_proba(100,:)>0)),'r')
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


%save exp_data spike_time50 spike_time100 perf_trial coherence RT
save datalast RT conf_level perf_trial perf mean_RT mean_RT_incor coherence neuron_traj antineuron_traj

beep
toc
