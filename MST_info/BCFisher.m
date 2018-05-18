function [FI varFI FInaive] = BCFisher(D1,D2,ds)
%% [FI varFI FInaive] = BCFisher(D1,D2,ds)
% Estimate Fisher Information (radians^-2) in a fine discrmination task.
%
% Inputs:
%
% D1 and D2 are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2. D1 and D2 must have the same size
% ds is in radians
%
% Outputs:
%
% FI, Fisher Information (radians^-2)
% varFI, variance of the estimator of Fisher Information
% FInaive, Fisher Information without bias correction (radians^-2)
%
% 
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

% These two are not good because prod(D1,2) could be Inf or NaN even if there's no Inf or NaN in D !!!! (HH20180509, test under Matlab2017)
% Also, the prod is too slow...
% ind1 = find(isnan(prod(D1,2)));  
% ind2 = find(isnan(prod(D2,2)));

ind1 = find(any(isnan(D1),2));  
ind2 = find(any(isnan(D2),2));  

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

Shat = squeeze(0.5*(nancov(D1)+nancov(D2))); % noise covariance
Fp = squeeze((nanmean(D1,1)-nanmean(D2,1))/ds); % signal
FInaive = (Fp*pinv(Shat)*Fp'); % FI, naive direct estimator
FI = FInaive*(2*T-2-N-1)/(2*T-2) - 2*N/(T*ds^2); % FI, bias corrected
varFI =  (2*FI^2)/(2*T-N-5) * (1 + 4*(2*T-3)/(T*FI*ds^2) + 4*N*(2*T-3)/(T*FI*ds^2)^2); % variance of FI bias corrected estimator

