

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization for lip.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



spikes_in_trial11 = zeros(N_in,nu_trial*2);
spikes_in_trial12 = zeros(N_in,nu_trial*2);

spikes_in21 = zeros(N_in,trial_dur);
spikes_in22 = zeros(N_in,trial_dur);

aux1_proba1 = zeros(N_out,trial_dur+2);
aux2_proba1 = zeros(N_out,trial_dur+2);


decision_ac = zeros(N_out,trial_dur+2);

aux1_proba2 = zeros(N_out,trial_dur+2);
aux2_proba2 = zeros(N_out,trial_dur+2);


spikes_in = cell(2,1);
proba_out = cell(2,1);
spikes_out = cell(2,1);
for m=1:2
   spikes_in{m} = zeros(N_out,trial_dur);
   proba_out{m} = zeros(N_out,trial_dur);
   spikes_out{m} = zeros(N_out,trial_dur);
end

proba_count = cell(2,nu_window);
spikes_count = cell(2,nu_window);
for m=1:2
for n=1:nu_window    
   proba_count{m,n} = zeros(N_out,nu_trial*2);
   spikes_count{m,n} = zeros(N_out,nu_trial*2);
   spikes_count_in{m,n} = zeros(N_in,nu_trial*2);
end
end

if shuffle_flag==1
   spikes_count_shuff = cell(2,nu_window);
   for m=1:2
   for n=1:nu_window    
      spikes_count_shuff{m,n} = zeros(N_out,nu_trial*2);
   end
   end
end

proba_out01  = zeros(N_out,pre_trial_dur);
proba_out02  = zeros(N_out,pre_trial_dur);
spikes_out01 = zeros(N_out,pre_trial_dur);
spikes_out02 = zeros(N_out,pre_trial_dur);

proba_out_av0__pretrial  = zeros(N_out,pre_trial_dur);
proba_out_av1  = zeros(N_out,trial_dur);

spikes_out_trial01 = zeros(N_out,nu_trial*2);
proba_out_trial01  = zeros(N_out,nu_trial*2);
spikes_out_trial02 = zeros(N_out,nu_trial*2);
proba_out_trial02  = zeros(N_out,nu_trial*2);

proba_out50_pre=zeros(nu_trial*2,pre_trial_dur);
proba_out1_pre=zeros(nu_trial*2,pre_trial_dur);
proba_out50=zeros(nu_trial*2,trial_dur);
proba_out1=zeros(nu_trial*2,trial_dur);

w_oi = zeros(N_out,N_in);
w_oi2 = zeros(N_out,N_in);
w_oo = zeros(N_out,N_out);

% << Network connections >>
for j=1:N_out
    w_oi(j,:) = g_w_oi/N_in *(exp(K_oin*(cos((pos_in-pos_out(j))/180 *2*pi)-1)));  % << MT input >>
    w_oi2(j,:) = g_w_oi2/N_in *(exp(K_oin2*(cos((pos_in-pos_out(j))/180 *2*pi)-1)));  % << Target input >>
    w_oo(j,:) = g_w_oo/N_out*...   % << LIP recurrent >>
                ((exp(K_oo*(cos((pos_out-pos_out(j))/180*2*pi)-1)))-...
                amp_i*(exp(K_oI*(cos((pos_out-pos_out(j))/180*2*pi)-1))))...
                + dc_w_oo;
end

% << MT correlation matrix >>
for j=1:N_in
    cov_mt(j,:) = var_mt*exp(K_cov_mt*(cos((pos_in-pos_in(j))/180 *2*pi)-1));
end
w_mt = real(sqrtm(cov_mt));

RT = zeros(nu_trial*2,1);


spikes_out50 = zeros(nu_trial*2,trial_dur);
spikes_out75 = zeros(nu_trial*2,trial_dur);
spikes_out100 = zeros(nu_trial*2,trial_dur);
last_spikes = zeros(N_out,2*nu_trial);
last_proba = zeros(N_out,2*nu_trial);
