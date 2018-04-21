function FI = Pcorrect2FI(PC,ds)
%% FI = Pcorrect2FI(PC,ds)
% Convert percent correct tto Fisher Information (radians^-2).
%
% Inputs:
%
% PC proportion correct
% ds is the difference between stimuli (radians)
%
% Outputs:
%
% FI is the fisher information
% 
% Copyright (c) 2015, Ingmar Kanitscheider and Ruben Coen Cagli. 
% All rights reserved.
% See the file LICENSE for licensing information.
%%

FI =  (2*norminv(PC)/ds)^2; 

