function [FIVAL, FITR] = RidgeRegLOOCV(D1,D2,ds,lambda)
%% [FIVAL, FITR] = RidgeRegLOOCV(D1,D2,ds,lambda)
% Computes Fisher information Using Ridge Regression (L2 regularization)
% with Leave One Out Cross Validation
%
% Inputs:
%
% D1 and D2 are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2. D1 and D2 must have the same size
% ds is in radians
% lambda is the weight on the penalty term
%
% Outputs:
% 
% FIVAL is the Fisher information on the Validation set (radians^-2)
% FITR is the Fisher information on the Training set (radians^-2)
%
%
% Ruben Coen-Cagli
% 2015-01

if(~exist('lambda'))
    lambda=0;
end
    
[ns1,nn]=size(D1);
[ns2,nn]=size(D2);

% NEED SAME # TRIALS IN EACH CONDITION
if ns1>ns2
    indtmp = randperm(ns1);
    D1 = D1(indtmp(1:ns2),:);
    ns1=ns2;
elseif ns2>ns1
    indtmp = randperm(ns2);
    D2 = D2(indtmp(1:ns1),:);
    ns2=ns1;
end
NR=ns1;
nsTR1 = NR-1;
nsTR2 = NR-1;
nsVAL1 = 1;
nsVAL2 = 1;

YTR = [-ones(nsTR1,1); ones(nsTR2,1)] * ds/2;

FITR = NaN(1,NR);
PRED1 = NaN(1,NR);
PRED2 = NaN(1,NR);
for k=1:NR

    indTR=[1:k-1 k+1:ns1];
    DTR1 = D1(indTR,:);
    DVAL1 = D1(k,:);
    DTR2 = D2(indTR,:);
    DVAL2 = D2(k,:);

    %  Remove the nanmean from Training set from both Training and Validation sets
    muTR = (nansum(DTR1)+nansum(DTR2))/(nsTR1+nsTR2);
    DTR1 = DTR1 - ones(nsTR1,1)*muTR;
    DTR2 = DTR2 - ones(nsTR2,1)*muTR;
    DVAL1 = DVAL1 - ones(nsVAL1,1)*muTR;
    DVAL2 = DVAL2 - ones(nsVAL2,1)*muTR;
    
    
    if lambda~=0
        w = ((nanmean(DTR2)-nanmean(DTR1)) / (0.5*(nancov(DTR1)+nancov(DTR2)) + lambda*eye(nn)) )';
    else % optimize lambda by cross-validation within the training set
        %%% fill this in, if you can wait two months for it to run through
    end

    
    % Estimate Fisher Information
    biasTR = (nanmean(DTR2)-nanmean(DTR1))*w/ds;
    varTR = w'*(nancov(DTR1)+nancov(DTR2))*w/2;
    FITR(k) = biasTR^2/(varTR*(nsTR1+nsTR2-2)/(nsTR1+nsTR2-4)) - 2/(0.5*(nsTR1+nsTR2)*ds^2); %*** 2014.04 RCC FI bias correction
    % Compute predictions 
    PRED1(k) = DVAL1*w;
    PRED2(k) = DVAL2*w;
end
FITR=nanmean(FITR);

% Compute Fisher Information from percent correct
pcVAL = (nanmean(PRED2 > 0)+nanmean(PRED1 < 0))/2;
FIVAL = Pcorrect2FI(pcVAL,ds);


