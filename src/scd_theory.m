%
%   File:      scd_theory.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/8

%   Description:
%       Obtains the theorical SCD of a single harmonic in additive noise.
%       xn = A*cos(2*pi*fc*n) + vn
%
%   Notes:
%       - Tested with Matlab 2007a


function [h f fa] = scd(A,fc,fs,N,res,var)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theorical SCD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = zeros(N,res);
f = [-N/2:N/2-1]*fs/N;
fa = [-res/2:res/2-1]*fs/res;
f1 = [-fliplr(fc); fc];     % frequency bins  
f2 = [-fliplr(2*fc); 2*fc]; % alpha bins (cyclic freq at -2fc and 2fc)
a1 = [fliplr(A) A];         % bins amplitude

for i = 1:length(f1)
    [m,ix] = min(abs(f-f1(i)));            % find index of the closest value to the frequency
    [m,i0] = min(abs(fa));                 % find index of the closest zero freq
    h(ix,i0) = (a1(i)^2)/4;                % theorical PSD (alpha = 0)
end

for i = 1:length(f2)
    [m,ix] = min(abs(fa-f2(i)));           % find index of the closest value to the frequency
    [m,i0] = min(abs(f));                  % find index of the closest zero freq
    h(i0,ix) = (a1(i)^2)/4;                % theorical PSD (alpha != 0)
end
