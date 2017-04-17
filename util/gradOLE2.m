%%%function foncOLE
%%% for amoeba usage

function [dE_dW]=gradOLE2(w)
global Y_targ X_train;

Y_train = w'* X_train;

%E_train = sum(sum( (Y_targ-Y_train).^2 ));


	dE_dW = -X_train*(Y_targ-Y_train)';
	%\dE_dW=dE_dW';


%grad_err = (sum(E_test(count-20:count-10))-sum(E_test(count-10:count)))/11;


