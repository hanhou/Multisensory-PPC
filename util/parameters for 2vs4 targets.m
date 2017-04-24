N_in = 100;
pos_in=[0:180/N_in:179.9999];
N_out = 100;
step_out = 180/N_out;
pos_out=[0:step_out:179.9999];
nu_targets = 4;
decis_thres = 56; %66 for 2 targets, 56 for 4
decis_flag = 1;

peak = 90;
pre_trial_dur= 50; %% must be 50 or more
trial_dur = 2000; %%all times are in ms

if pre_trial_dur>trial_dur
   error('pre_trial_dur cannot be longer than trial_dur'); 
end

delta_t = 1e-3; %size of time bin in second
% difference in s between the two test stimuli
step_test =2;


% Parameters for MT from Mazurek and Shadlen, except width. K_cov_mt  
% and var_mt controls the correlations. 
r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;
K_cov_mt = 2;
var_mt = 1e-5;

% parameter of the visually evoked activity
max_rate_in2 = 42;
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
g_w_oi = 70;

% Visual evoked weights
g_w_oi2= 40;
%%% drop in attention to visual target once motion stimulus appears.
att_gain2 = 0;


g_w_oo = 21;
amp_i = 0.0; 
dc_w_oo = -0.12;

ext_noise=0;

shuffle_flag = 0;
info_out_flag = 1;
info_in_flag = 1;

%%% duration of the window used to compute Fisher
wind_size = 50;
wind_offset = 0;
nu_window = 4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialization

