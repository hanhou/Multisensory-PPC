function [FIVAL, FITR, FIVAL_E] = EarlyStoppingShufDiag(D1,D2,ds,fracTR,fracTE,NR,AVG)
%% [FIVAL, FITR, FIVAL_E] = EarlyStoppingShufDiag(D1,D2,ds,fracTR,fracTE,NR,AVG)
% Computes Fisher information Using Early Stopping
%*** I_shuffle: Train on shuffled data , validate on shuffled data ***%
%*** I_diag: Train on shuffled data , validate on original data ***%
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
% FIVAL is I_shuffle on the Validation set (radians^-2)
% FITR is I_shuffle on the Training set (radians^-2)
% FIVAL_E is I_diag on the Validation set (radians^-2)
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
    AVG=0
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

nsE1 = nsVAL1;
nsE2 = nsVAL2;
nsEVAL1 = nsVAL1;
nsEVAL2 = nsVAL2;

FIVAL_E = NaN(1,NR);
FIVAL = NaN(1,NR);
FITR = NaN(1,NR);
for k=1:NR
    
    idx=randperm(ns1);
    DTR1 = D1(idx(1:nsTR1),:);
    DTE1 = D1(idx(nsTR1+1:nsTR1+nsTE1),:);
    DVAL1 = D1(idx(nsTR1+nsTE1+1:ns1),:);
    
    idx=randperm(ns2);
    DTR2 = D2(idx(1:nsTR2),:);
    DTE2 = D2(idx(nsTR2+1:nsTR2+nsTE2),:);
    DVAL2 = D2(idx(nsTR2+nsTE2+1:ns2),:);
    
    EVAL1 = DVAL1;
    EVAL2 = DVAL2;
    
    % SHUFFLE TRAINING, TEST AND VALIDATION SETS
    if 1%(nsTR1<nn || nsTE1<nn || nsVAL1<nn ||nsTR2<nn || nsTE2<nn || nsVAL2<nn)
        for n=1:nn
            DTR1(:,n) = DTR1(randperm(nsTR1),n);
            DTR2(:,n) = DTR2(randperm(nsTR2),n);
            DTE1(:,n) = DTE1(randperm(nsTE1),n);
            DTE2(:,n) = DTE2(randperm(nsTE2),n);
            DVAL1(:,n) = DVAL1(randperm(nsVAL1),n);
            DVAL2(:,n) = DVAL2(randperm(nsVAL2),n);
        end
    else %***** Use this if you want to make sure each neuron uses a different trial -- need more trials than neurons! *****%
        tmpDTE1=DTE1;
        tmpDTE2=DTE2;
        tmpDTR1=DTR1;
        tmpDTR2=DTR2;
        tmpDVAL1=DVAL1;
        tmpDVAL2=DVAL2;
        for t=1:nsTR1
            indTR = randperm(nsTR1);
            for n=1:nn
                tmpDTR1(t,n) = DTR1(indTR(n),n);
            end
        end
        for t=1:nsTR2
            indTR = randperm(nsTR2);
            for n=1:nn
                tmpDTR2(t,n) = DTR2(indTR(n),n);
            end
        end
        for t=1:nsTE1
            indTE = randperm(nsTE1);
            for n=1:nn
                tmpDTE1(t,n) = DTE1(indTE(n),n);
            end
        end
        for t=1:nsTE2
            indTE = randperm(nsTE2);
            for n=1:nn
                tmpDTE2(t,n) = DTE2(indTE(n),n);
            end
        end
        for t=1:nsVAL1
            indVAL = randperm(nsVAL1);
            for n=1:nn
                tmpDVAL1(t,n) = DVAL1(indVAL(n),n);
            end
        end
        for t=1:nsVAL2
            indVAL = randperm(nsVAL2);
            for n=1:nn
                tmpDVAL2(t,n) = DVAL2(indVAL(n),n);
            end
        end
        DTE1=tmpDTE1;
        DTE2=tmpDTE2;
        DTR1=tmpDTR1;
        DTR2=tmpDTR2;
        DVAL1=tmpDVAL1;
        DVAL2=tmpDVAL2;
    end
    
    %  Remove the nanmean from Training set from both Training and Test sets
    muTR = (nansum(DTR1)+nansum(DTR2))/(nsTR1+nsTR2);
    DTR1 = DTR1 - ones(nsTR1,1)*muTR;
    DTR2 = DTR2 - ones(nsTR2,1)*muTR;
    DTE1 = DTE1 - ones(nsTE1,1)*muTR;
    DTE2 = DTE2 - ones(nsTE2,1)*muTR;
    DVAL1 = DVAL1 - ones(nsVAL1,1)*muTR;
    DVAL2 = DVAL2 - ones(nsVAL2,1)*muTR;
    
    EVAL1 = EVAL1 - ones(nsEVAL1,1)*muTR;
    EVAL2 = EVAL2 - ones(nsEVAL2,1)*muTR;
    
    % Compute residuals for the training sets
    sbarTR = ds/2*(nsTR2-nsTR1)/(nsTR2+nsTR1);
    mupTR = (nansum(DTR2)-nansum(DTR1))/(nsTR2+nsTR1)*ds/2;
    
    muTE = (nansum(DTE1)+nansum(DTE2))/(nsTE1+nsTE2);
    mupTE = (nansum(DTE2)-nansum(DTE1))/(nsTE2+nsTE1)*ds/2;
    
    COVTR=(DTR1'*DTR1+DTR2'*DTR2)/(nsTR1+nsTR2);
    M2TE = (DTE1'*DTE1+DTE2'*DTE2)/(nsTR1+nsTR2);
    
    
    
    %********** Optimization loop
    dt=1/10/nanmax(eig(COVTR));
    dETEdt = -Inf;
    iters(k)=0;
    w=zeros(nn,1);  %vector of weights
    dETRdw = COVTR*w - mupTR';
    while(dETEdt<0 & iters(k) < maxiters)
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
    biasTR = (nanmean(DTR2)-nanmean(DTR1))*w/ds;
    varTR = w'*(nancov(DTR1)+nancov(DTR2))*w/2;
    FITR(k) = biasTR^2/(varTR*(nsTR1+nsTR2-2)/(nsTR1+nsTR2-4)) - 2/(0.5*(nsTR1+nsTR2)*ds^2); %*** 2014.04 RCC FI bias correction
    
    biasVAL = (nanmean(DVAL2)-nanmean(DVAL1))*w/ds;
    varVAL = w'*(nancov(DVAL1)+nancov(DVAL2))*w/2;
    FIVAL(k) = biasVAL^2/(varVAL*(nsVAL1+nsVAL2-2)/(nsVAL1+nsVAL2-4)) - 2/(0.5*(nsVAL1+nsVAL2)*ds^2); %*** 2014.04 RCC FI bias correction
    
    biasVAL_E = (nanmean(EVAL2)-nanmean(EVAL1))*w/ds;
    varVAL_E = w'*(nancov(EVAL1)+nancov(EVAL2))*w/2;
    FIVAL_E(k) = biasVAL_E^2/(varVAL_E*(nsEVAL1+nsEVAL2-2)/(nsEVAL1+nsEVAL2-4)) - 2/(0.5*(nsEVAL1+nsEVAL2)*ds^2); %*** 2014.04 RCC FI bias correction
end

if AVG % average results across cross-validation splits
    FIVAL=nanmean(FIVAL);
    FITR=nanmean(FITR);
    FIVAL_E = nanmean(FIVAL_E);
end

%%
