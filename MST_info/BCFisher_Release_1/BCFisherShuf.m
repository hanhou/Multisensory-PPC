function [FI FInaive] = BCFisherShuf(D1,D2,ds)
%% [FI FInaive] = BCFisherShuf(D1,D2,ds)
% Estimate Fisher Information (radians^-2) in a fine discrmination task,
% assuming an independent population. 
% Known as I_shuffle.
%
% Inputs:
%
% D1 and D2 are data matrices (trials x neurons) under the two
% conditions s=-ds/2, s=+ds/2.
% ds is in radians
%
% Outputs:
%
% FI, Fisher Information (radians^-2)
% FInaive, Fisher Information without bias correction (radians^-2)
%
%*** NOTE ***%
% This analytical correction is more accurate (and faster) than shuffling
% the data and then using BCFisher(D1shuf, D2shuf, ds).
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

V1 = squeeze(0.5*(nanvar(D1)+nanvar(D2))); % variance
Fp1 = squeeze((nanmean(D1)-nanmean(D2)).^2./(ds^2)); % signal
FI1= (Fp1./V1); % single-neuron FI, naive estimator

ind = isfinite(FI1); % remove neurons with 0 variance (eg due to no spikes)
N = sum(ind);
FInaive = nansum(FI1(ind)); % population FI, naive estimator
FI = FInaive*(T-2)/(T-1) - 2*N/(T*ds^2); % FI, bias corrected

%%