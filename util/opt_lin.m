%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plain vanilla gradient descent for a quadratic cost function 
% using early stopping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W_out = opt_lin(learn_rate,Y_targ,X_train,X_test,W);

max_trial = 300000;


E_train = zeros(1,max_trial);
E_test = zeros(1,max_trial);

Y_test = W'*X_test;
Y_test_ini = Y_test;
E_test(1:20) = sum(sum( (Y_targ-Y_test).^2 ));
grad_err =0;

count =21;

while count<max_trial+1
   Y_train = W'* X_train;
   Y_test = W'* X_test;

   E_train(count) = sum(sum( (Y_targ-Y_train).^2 ));
   E_test(count) =  sum(sum( (Y_targ-Y_test).^2 ));
   if mod(count,100)==0
      fprintf('count:  %d E_train: %e    E_test: %e \n',count,E_train(count),E_test(count));
      save w_ole_temp W;    
   end
   stop_count = count;
   grad_err = (sum(E_test(count-20:count-10))-sum(E_test(count-10:count)))/11;
   if (grad_err<0)
       fprintf('Early stopping: %d\n',count);
       fprintf('Difference: %e\n',grad_err);
       count =max_trial;
   else 
      dE_dW = -X_train*(Y_targ-Y_train)';

      W = W - learn_rate* dE_dW;
   end
  count = count+1;
end

if stop_count==max_trial+1
    fprintf('Maximum number of trials reached: %d \n',stop_count-1);
end

W_out = W;

sampled = E_train(1,1:100:max_trial);
save error sampled;


subplot(221)
plot(Y_train);
hold on
plot(Y_test,'c');
plot(Y_targ,'r');
plot(Y_test_ini,'g');
hold off

subplot(222)
plot(log10(E_train(11:stop_count)),'r');
hold on;
plot(log10(E_test(11:stop_count)));
hold off;

subplot(223);
plot(W_out,'r');
