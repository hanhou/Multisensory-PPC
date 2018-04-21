function [FIE VARE BIASE] = BCFisherCross(D1,D2,E1,E2,ds)
%% [FIE VARE BIASE varFIE] = BCFisherCross(D1,D2,E1,E2,ds)
% Estimate Fisher Information (radians^-2) in a fine discrmination task,
% assuming suboptimal decoder, i.e. optimized to response
% statistics of population D and tested on population E.
%
% Inputs: 
%
% D1 (E1) and D2 (E2) are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2.
% ds is in radians
%  
% Outputs:
%
% FIE, Fisher Information (radians^-2)
% VARE, decoder variance on test data: wD*SigmaE*wD'
% BIASE, decoder bias derivative on test data: (wD*muprimeE)
%
% Copyright (c) 2015, Ingmar Kanitscheider and Ruben Coen Cagli. 
% All rights reserved.
% See the file LICENSE for licensing information.
% 
% For derivation, see:
% Kanitscheider*, Coen-Cagli*, Kohn, Pouget. "Measuring Fisher information 
% accurately in correlataed neural populations". PLoS Comp Biol
% 2015

%% Preliminary checks

[T N] = size(D1);% T=trials; N=neurons

if(size(D2,1)~=T)
    warning('The two datasets have different numbers of trials. Results may be inaccurate.')
end
if(size(D2,2)~=N)
    error('The two datasets have different numbers of neurons.')
end

% If there are NaNs, subsample to make D1 and D2 same size
ind1 = find(isnan(prod(D1,2)));
ind2 = find(isnan(prod(D2,2)));
N1=numel(ind1);
N2=numel(ind2);
if(N1>0 || N2>0)
    warning('Removing rows with NaNs. Subsampling to make the two datasets equal size')
    T = min([N1 N2]);
    indT1 = randperm(numel(ind1),T);
    indT2 = randperm(numel(ind2),T);
    D1 = D1(indT1,:);
    D2 = D2(indT2,:);
end

%% Compute outputs

% sample estimates of mean, derivative, covariance
emuD1 = nanmean(D1,1);
emuD2 = nanmean(D2,1);
emuE1 = nanmean(E1,1);
emuE2 = nanmean(E2,1);
emupD = (emuD1-emuD2)/ds; % tuning curve derivatives
emupE = (emuE1-emuE2)/ds;
eCD = 0.5*(nancov(D1)+nancov(D2)); % sample covariances
eCE = 0.5*(nancov(E1)+nancov(E2));
einvCD = pinv(eCD); % naive inverse covariance
invCD = einvCD*(2*T-2-N-1)/(2*T-2); % bias-corrected inverse covariance

% bias-corrected estimate of information in D (optimal decoder)
FID = (emupD*invCD*emupD') - 2*N/(T*ds^2); % FI, bias corrected
   
% bias-corrected estimate of crossed information in E (suboptimal decoder)
d = (2*T-N-2)*(2*T-N-5);
tr = trace(eCE*invCD); % bias-corrected, invCD
lambda = (2/(T*ds^2)) * (tr *(1+((2*T-N-1)+N*(2*T-N-3))/d));
rho = FID * tr * (2*T-N-3)/d; % bias-corrected, FID    
VARnum = ((emupD*invCD*eCE*invCD*emupD') - lambda - rho) / ((1+(2*T-N-1)/d)); % bias-corrected
VARden = FID^2;
VARE = VARnum / VARden;
BIASE = (emupD*invCD*emupE') / FID;
FIE = (emupD*invCD*emupE')^2 / VARnum; % BIASE^2/VARE;

%%