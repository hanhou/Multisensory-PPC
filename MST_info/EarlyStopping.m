function [FIVAL, FITR] = EarlyStopping(D1,D2,ds,fracTR,fracTE,NR,AVG)
%% [FIVAL, FITR] = EarlyStopping(D1,D2,ds,fracTR,fracTE,NR,AVG)
% Computes Fisher information Using Early Stopping
%
% Inputs:
%
% D1 and D2 are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2. D1 and D2 must have the same size
% ds is in radians
% fracTR is the fraction of the data to be used for Training (typically 1/3)
% fracTE is the fraction of the data to be used for Stopping (typically 1/3)
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
% Code modified after Jeff Beck's early stopping 
% explained in SI of (Moreno-Bote et al., Nat Neurosci, 2014)
%

if(~exist('AVG'))
    AVG=0;
end
maxiters=10000;

[ns1,nn]=size(D1);
[ns2,nn]=size(D2);

% Generate Training and Testing sets
nsTR1=floor(ns1*fracTR);
nsTR2=floor(ns2*fracTR);
nsTE1=floor(ns1*fracTE);
nsTE2=floor(ns1*fracTE);

nsVAL1=ns1-nsTE1-nsTR1;
nsVAL2=ns2-nsTE1-nsTR2;

% Initialize output variables
FIVAL = NaN(1,NR);
FITR = NaN(1,NR);
parfor k=1:NR

    idx=randperm(ns1);
    DTR1 = D1(idx(1:nsTR1),:);
    DTE1 = D1(idx(nsTR1+1:nsTR1+nsTE1),:);
    DVAL1 = D1(idx(nsTR1+nsTE1+1:ns1),:);
    
    idx=randperm(ns2);
    DTR2 = D2(idx(1:nsTR2),:);
    DTE2 = D2(idx(nsTR2+1:nsTR2+nsTE2),:);
    DVAL2 = D2(idx(nsTR2+nsTE2+1:ns2),:);

    %  Remove the mean of Training set from both Training and Test sets
    muTR = (sum(DTR1)+sum(DTR2))/(nsTR1+nsTR2);
    DTR1 = DTR1 - ones(nsTR1,1)*muTR;
    DTR2 = DTR2 - ones(nsTR2,1)*muTR;
    DTE1 = DTE1 - ones(nsTE1,1)*muTR;
    DTE2 = DTE2 - ones(nsTE2,1)*muTR;
    DVAL1 = DVAL1 - ones(nsVAL1,1)*muTR;
    DVAL2 = DVAL2 - ones(nsVAL2,1)*muTR;

    % Compute residuals for the training sets
    sbarTR = ds/2*(nsTR2-nsTR1)/(nsTR2+nsTR1);
    mupTR = (sum(DTR2)-sum(DTR1))/(nsTR2+nsTR1)*ds/2;

    muTE = (sum(DTE1)+sum(DTE2))/(nsTE1+nsTE2);
    mupTE = (sum(DTE2)-sum(DTE1))/(nsTE2+nsTE1)*ds/2;

    COVTR=(DTR1'*DTR1+DTR2'*DTR2)/(nsTR1+nsTR2); % DTR1/2 are residuals
    M2TE = (DTE1'*DTE1+DTE2'*DTE2)/(nsTR1+nsTR2); % Should be nsTE1/2???

    
    %********** Optimization loop <Verified by HH.>
    dt=1/10/max(eig(COVTR)); % define stepsize
    dETEdt = -Inf;
    iters(k)=0;
    w=zeros(nn,1);  % initial vector of weights
    dETRdw = COVTR*w - mupTR';
    while(dETEdt<0 && iters(k) < maxiters)
       iters(k)=iters(k)+1;
       w = w - dt*dETRdw; % update w

       dETRdw = COVTR*w - mupTR';
       dETEdw = M2TE*w - mupTE' + sbarTR*muTE';  

       dETEdt = -dETRdw'*dETEdw;

    end
    if(iters(k)==maxiters)
       fprintf(strcat('Max iters reached -- run ',num2str(k),'\n'))
    end
    %********** End of optimization loop
    
    
    % Estimate Fisher Information     
    biasTR = (mean(DTR2)-mean(DTR1))*w/ds; 
    varTR = w'*(cov(DTR1)+cov(DTR2))*w/2;
% % %     FITR(k) = biasTR^2/varTR; %*** Original version by Jeff Beck
    FITR(k) = biasTR^2/(varTR*(nsTR1+nsTR2-2)/(nsTR1+nsTR2-4)) - 2/(0.5*(nsTR1+nsTR2)*ds^2); 
    %*** 2014.04 RCC added: bias correction (Where is this last term from??? HH)
    
    biasVAL = (mean(DVAL2)-mean(DVAL1))*w/ds; 
    varVAL = w'*(cov(DVAL1)+cov(DVAL2))*w/2;
% % %     FIVAL(k) = biasVAL^2/varVAL; %*** Original version by Jeff Beck
    FIVAL(k) = biasVAL^2/(varVAL*(nsVAL1+nsVAL2-2)/(nsVAL1+nsVAL2-4)) - 2/(0.5*(nsVAL1+nsVAL2)*ds^2); %*** 2014.04 RCC added: bias correction   
    
end


if AVG % average results across cross-validation splits
    FIVAL=nanmean(FIVAL);
    FITR=nanmean(FITR);
end
