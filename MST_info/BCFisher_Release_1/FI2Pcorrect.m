function PC = FI2Pcorrect(FI,ds)
%% PC = FI2Pcorrect(FI,ds)
% Convert Fisher Information (radians^-2) to percent correct.
%
% Inputs:
%
% FI is the fisher information
% ds is the difference between stimuli (radians)
%
% Outputs:
%
% PC proportion correct
% 
% Copyright (c) 2015, Ingmar Kanitscheider and Ruben Coen Cagli. 
% All rights reserved.
% See the file LICENSE for licensing information.
%%

STDest = 1./(FI*(pi/180)^2).^.5; % STD of optimal unbiased estimator, in degrees
PC = normcdf(0.5*ds./STDest);

