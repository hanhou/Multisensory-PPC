%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program loads real data and uses gradient descent
% to find the weights of the OLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
randn('state',sum(100*clock))

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
global Y_targ X_train;
global nu_unit
global target
global nu_trial
K_f=1.5;
K_cov=3;
spont=1;
max_mean=39;
b=0;
global ID
sig_prior=0.000000000005;
global cov_act
global der_act
global theta_pref
%%% learning rate. Lower value used for input patterns
%learn_rate = 6e-9
%learn_rate = 4e-12
load ..\knill\w_knill
%%%% init_weights=0 means that initial weights are picked at random
%%%% otherwise they are set to the weights loaded in the previous line
init_weights=1;
%%% conjugate gradient flag
cg_flag=1;
%%% gradient descent flag
gd_flag=0;
%%% Shuffled test: used shuffle data for training
shuffled=0;
%%% Swap data: exchange training data for test data
swap_data=0;
%%% subsample flag. Set to 1 is only a subset of the neurons are used.
subsample=0;
%%% Subsample some of the data. Must be 1 or a multiple of 4.
neuron_step=1;
%%% total number of sample. set to maximum if equal to 0.
total_subsample=0;
%%% first subsampled neurons
first_neuron = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../knill/knilldata
%load peggy\sharpen1\datain895905;
%load peggy\sharpen1\datain895905_shuffled;
%data1=data895;
%data2=data905;

%%% step size: half the angle between the two orientations
step_test=2;

[nu_trial nu_unit] = size(data1);
if subsample==1 
   if total_subsample==0
       last_neuron = nu_unit;
   else
      last_neuron = first_neuron+(total_subsample-1)*neuron_step;
   end
   if last_neuron>nu_unit
       error('Total subsample is too high\n');
   end
   aux1=data1;
   aux2=data2;
   clear data1 data2 
   data1 = aux1(:,first_neuron:neuron_step:last_neuron); 
   data2 = aux2(:,first_neuron:neuron_step:last_neuron); 
   
   aux3 = act_train_shuffled;
   aux4 = act_test_shuffled;
   clear act_train_shuffled act_test_shuffled
   act_train_shuffled=aux3(:,first_neuron:neuron_step:last_neuron); 
   act_test_shuffled=aux4(:,first_neuron:neuron_step:last_neuron); 

   
   [nu_trial nu_unit] = size(data1);
   
   fprintf('\nCAUTION: sampling every %d neurons\n',neuron_step); 
   fprintf('\nCAUTION: using only %d neurons\n',nu_unit); 
   fprintf('\nCAUTION: using neurons %d to %d\n\n',first_neuron,last_neuron); 
   clear aux1 aux2 aux3 aux4;
end 


%% next 2 line: compute the preferred directions
step=180/(nu_unit/4);
for k=1:nu_unit/4
    for l=1:4
       theta_pref((k-1)*4+l)= (-90+step*(k-1))/180*pi;
    end
end

nu_trial=nu_trial/2;
% aux_nu_trial alows to select the first aux_nu_trial of the data set
aux_nu_trial = nu_trial ;
%aux_nu_trial = 500;

if(aux_nu_trial~=nu_trial)
    fprintf('\nCAUTION: using only a subset of all trials\n\n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% reorganizes data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

act_train(1:aux_nu_trial,:)= data1(1:aux_nu_trial,:); 
act_train(aux_nu_trial+1:2*aux_nu_trial,:) = data2(1:aux_nu_trial,:);
act_test(1:aux_nu_trial,:)= data1(nu_trial+1:nu_trial+aux_nu_trial,:); 
act_test(aux_nu_trial+1:2*aux_nu_trial,:) = data2(nu_trial+1:nu_trial+aux_nu_trial,:);
%clear data1 data2;

if shuffled==1
   act_train(1:aux_nu_trial,:)= data1_shuffled(1:aux_nu_trial,:); 
   act_train(aux_nu_trial+1:2*aux_nu_trial,:) = data2_shuffled(1:aux_nu_trial,:);
   act_test(1:aux_nu_trial,:)= data1_shuffled(nu_trial+1:nu_trial+aux_nu_trial,:); 
   act_test(aux_nu_trial+1:2*aux_nu_trial,:) = data2_shuffled(nu_trial+1:nu_trial+aux_nu_trial,:);
   fprintf('\nCAUTION: using shuffled data\n\n');
end

if swap_data==1
   junk = act_test;
   act_test =act_train;
   act_train = junk;
   fprintf('\nCAUTION:Swapping data\n');
   clear junk;
end

Y_targ(1:aux_nu_trial,1) = step_test; 
Y_targ(aux_nu_trial+1:2*aux_nu_trial,1) = -step_test;

nu_trial=aux_nu_trial;
fprintf('Number of trials: %d\n',nu_trial);
display('Data loaded');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% compute mean and derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_act = mean(act_train(1:nu_trial,:));
cov_act = cov(act_train(1:nu_trial,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Classify with K nearest neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%K=1;
%fprintf('Starting K Nearest neighbor. K=%d\n\n',K);
%labels =  KNearestNeighbors(K,act_test.^2,act_train.^2,Y_targ);
%perc_right_test_KNB = sum(max(Y_targ)*((labels>0)*2-1)==Y_targ')/(nu_trial*2)*100
%fprintf('Done with K Nearest neighbor\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find optimum weights using grad descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if init_weights==0
   fprintf('Initializing weights\n');
   w_ole_ini = rand(nu_unit+1,1)*1e-8; 
   w_ole_ini = w_ole_ini-mean(w_ole_ini);
   w_ole_ini(nu_unit+1)=0; 
else
   fprintf('Using preloaded weights\n');
   [size_w_ole junk2] = size(w_ole);
   if size_w_ole == nu_unit
      w_ole(nu_unit+1,1)=0;
   end
   w_ole_ini = w_ole;
end

act_train(:,nu_unit+1)=ones(nu_trial*2,1);
act_test(:,nu_unit+1)=ones(nu_trial*2,1);
act_train_zero_mean = act_train-ones(nu_trial*2,1)*mean(act_train);
clear act_train;
act_test_zero_mean = act_test-ones(nu_trial*2,1)*mean(act_test);
clear act_test;

if gd_flag==1
   disp('Starting batch linear regression');
      w_ole = opt_lin(learn_rate,Y_targ',act_train_zero_mean',act_test_zero_mean',w_ole_ini);
      save w w_ole
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find optimum weights with conjugate gradient descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cg_flag==1
   disp('Minimization using conjugate gradients methods');
    X_train=act_train_zero_mean';
    Y_targ=Y_targ';
    gtol = eps; 
    ftol=1e-14;
	iters=1;
	c1=0.001; c2=0.6; %%necessary 0<c1<c2<1 Wolfe line search
    info_test_min=0;
    count=0;
    k=0;
    max_iter=10000;
    while k<max_iter
	   [w_ole_cg,err,foo1,foo2,foo3,exitflag] = ...
           cg_pr('foncgradOLE2','gradOLE2',w_ole_ini,gtol,ftol,c1,c2,iters,0); % Polak-Ribiere conjugate gradient method (without restart) 
       info_test = fisher(w_ole_cg,act_test_zero_mean,step_test);
       fprintf('Iter: %d  I_test: %6.4f  I_test max: %6.4f  count:%d\n',...
           k*iters,info_test,info_test_min,count);
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
    end
    Y_targ=Y_targ';
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% testing the weights on data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info_test = fisher(w_ole,act_test_zero_mean,step_test);
info_train = fisher(w_ole,act_train_zero_mean,step_test);

theta_hat_test =  w_ole'* act_test_zero_mean';
theta_hat_train =  w_ole'* act_train_zero_mean';

perc_right_train = sum(max(Y_targ)*((theta_hat_train>0)*2-1)==Y_targ')/(nu_trial*2)*100;
perc_right_test= sum(max(Y_targ)*((theta_hat_test>0)*2-1)==Y_targ')/(nu_trial*2)*100;

pred_d_prime = erfinv(2*(perc_right_test/100-0.5))*2^0.5;
pred_info= ((pred_d_prime*2)/step_test)^2;
d_prime= (info_test^0.5*step_test)/2;

std_test= ((var(theta_hat_test(1:nu_trial))+var(theta_hat_test(nu_trial+1:nu_trial*2)))/2)^.5;
pred_perc = ( (0.5+(erf(d_prime/2^.5)/2)) )*100;
pred_thre = 2*(1.35/info_test^.5);

fprintf('\n');
fprintf('Info train          = %6.4f  %6.4f\n',info_train^.5,info_train);
fprintf('Info test           = %6.4f  %6.4f\n\n',info_test^.5,info_test);
fprintf('Predicted info Test = %6.4f  %6.4f\n\n',pred_info^.5,pred_info);

fprintf('Percentage correct on training data : %3.2f\n',perc_right_train);
fprintf('Percentage correct on test data     : %3.2f\n',perc_right_test);
fprintf('Predicted Perc correct on test data : %3.2f\n',pred_perc);
fprintf('Predicted discrimination threshold  : %3.2f deg\n\n',pred_thre);
 
fprintf('d_prime      from info        = %6.4f\n',d_prime);
fprintf('pred d_prime from performance = %6.4f\n\n',pred_d_prime);

aux_cov=diag(cov_act).^.5*(diag(cov_act).^.5)';
cor_coeff=cov_act./(aux_cov+1e-15);
cor_coeff=(cor_coeff<0.99999).*cor_coeff;
fprintf('Range of cor coeff = %6.4f   %6.4f   Mean: %6.4f   Var:%6.4f\n\n',...
    min(min(cor_coeff)), max(max(cor_coeff)),mean(mean(cor_coeff)),mean(var(cor_coeff))^.5);


test=zeros(nu_unit,nu_unit);
row=[1:nu_unit]'*ones(1,nu_unit);
column=row';
test=(abs(row-column)<nu_unit/16) | (abs(row-column)>15*nu_unit/16);
aux_coeff1=test.*cor_coeff;
aux_coeff2=(1-test).*cor_coeff;


disp('Preferred orientations less than 11 degrees apart');
fprintf('Range of cor coeff = %6.4f   %6.4f   Mean: %6.4f   Std:%6.4f\n\n',...
    min(min(aux_coeff1)), max(max(aux_coeff1)),sum(sum(aux_coeff1))/sum(sum(test))...
    ,mean(var(aux_coeff1))^.5);
disp('Preferred orientations more than 11 degrees apart');
fprintf('Range of cor coeff = %6.4f   %6.4f   Mean: %6.4f   Std:%6.4f\n\n',...
    min(min(aux_coeff2)), max(max(aux_coeff2)),sum(sum(aux_coeff2))/sum(sum(1-test))...
    ,mean(var(aux_coeff2))^.5);


if shuffled==1
   clear act_train act_test;
   act_train(1:aux_nu_trial,:)= data1(1:aux_nu_trial,:); 
   act_train(aux_nu_trial+1:2*aux_nu_trial,:) = data2(1:aux_nu_trial,:);
   act_test(1:aux_nu_trial,:)= data1(nu_trial+1:nu_trial+aux_nu_trial,:); 
   act_test(aux_nu_trial+1:2*aux_nu_trial,:) = data2(nu_trial+1:nu_trial+aux_nu_trial,:); 
   
   act_train(:,nu_unit+1)=ones(nu_trial*2,1);
   act_test(:,nu_unit+1)=ones(nu_trial*2,1);
   act_train_zero_mean = act_train-ones(nu_trial*2,1)*mean(act_train);
   act_test_zero_mean = act_test-ones(nu_trial*2,1)*mean(act_test);

   info_test = fisher(w_ole,act_test_zero_mean,step_test);
   info_train = fisher(w_ole,act_train_zero_mean,step_test);

   fprintf('Info train diag     = %8.2f  %6.4f\n',info_train^.5,info_train);
   fprintf('Info test diag      = %8.2f  %6.4f\n',info_test^.5,info_test);  
end
    
    
%%%%%%%%%%%%%%%%%%%%%%
%%%% plots
%%%%%%%%%%%%%%%%%%%%%%
subplot(221)
%surf(cov_act(1:20:nu_unit,1:20:nu_unit));
%surf(cor_coeff(1:10:nu_unit,1:10:nu_unit));
%title('Estimated covariance')
%axis([1 51 1 51 -25 50])
%caxis([-25 50])
%view(-60,44);

subplot(222)
plot(mean_act,'r');
axis tight;

subplot(223)
plot(Y_targ)
hold on
plot(theta_hat_test,'g')
hold off

subplot(224)
plot(w_ole,'g');
title('OLE weights')
axis tight;


