function [FIVAL, FITR] = VBLogReg(D1,D2,ds,fracTR,NR,AVG)
%% [FIVAL, FITR] = VBLogReg(D1,D2,ds,fracTR,NR,AVG)
% Computes Fisher information Using Variational Bayes Logistic Regression
% with automatic relevance detection
%
% Inputs:
%
% D1 and D2 are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2. D1 and D2 must have the same size
% ds is in radians
% fracTR is the fraction of the data to be used for Training
% NR is number of Runs with different permutations of the data
% AVG set to 'true' if results should be averaged across cross-validation splits
%
% Outputs:
% 
% FIVAL is the Fisher information on the Validation set (radians^-2)
% FITR is the Fisher information on the Training set (radians^-2)
%
%
% Ruben Coen-Cagli
% 2014-07
% 2015-01
%
% Variational Bayes Logistic Regression code provided by Jan Drugowitsch
% (Drugowitsch, arxiv:1310.5438 [stat.ML], 2014)
%


if(~exist('AVG'))
    AVG=0;
end
    
[ns1,nn]=size(D1);
[ns2,nn]=size(D2);

% Generate Training and Validation sets
nsTR1=floor(ns1*fracTR);
nsTR2=floor(ns2*fracTR);
nsVAL1=ns1-nsTR1;
nsVAL2=ns2-nsTR2;

YTR = [-ones(nsTR1,1); ones(nsTR2,1)];

FIVAL = NaN(1,NR);
FITR = NaN(1,NR);
for k=1:NR
    idx=randperm(ns1);
    DTR1 = D1(idx(1:nsTR1),:);
    DVAL1 = D1(idx(nsTR1+1:ns1),:);
    
    idx=randperm(ns2);
    DTR2 = D2(idx(1:nsTR2),:);
    DVAL2 = D2(idx(nsTR2+1:ns2),:);

    %  Remove the nanmean from Training set from both Training and Test sets
    muTR = (nansum(DTR1)+nansum(DTR2))/(nsTR1+nsTR2);
    DTR1 = DTR1 - ones(nsTR1,1)*muTR;
    DTR2 = DTR2 - ones(nsTR2,1)*muTR;
    DVAL1 = DVAL1 - ones(nsVAL1,1)*muTR;
    DVAL2 = DVAL2 - ones(nsVAL2,1)*muTR;
    
    [w, ~,~,~,~,~] = vb_logit_fit_ard([DTR1;DTR2], YTR);
    
    % Estimate Fisher Information
    biasTR = (nanmean(DTR2)-nanmean(DTR1))*w/ds;
    varTR = w'*(nancov(DTR1)+nancov(DTR2))*w/2;
    FITR(k) = biasTR^2/(varTR*(nsTR1+nsTR2-2)/(nsTR1+nsTR2-4)) - 2/(0.5*(nsTR1+nsTR2)*ds^2); %*** 2014.04 RCC FI bias correction

    tmpbias = (nanmean(DVAL2)-nanmean(DVAL1))*w/ds;
    tmpvar = w'*(nancov(DVAL1)+nancov(DVAL2))*w/2;
    FIVAL(k) = tmpbias^2/(tmpvar*(nsVAL1+nsVAL2-2)/(nsVAL1+nsVAL2-4)) - 2/(0.5*(nsVAL1+nsVAL2)*ds^2); %*** 2014.04 RCC FI bias correction
end


if AVG 
    FIVAL=nanmean(FIVAL);
    FITR=nanmean(FITR);
end

