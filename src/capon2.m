%
%   File:      capon2.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17
%
%   Description:
%       Estimates the cross PSD of two signals using the Capon2 (NMLM) estimator.
%       This implementation obtains the PSD by iterating over the range of frequencies of interes.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9

function [C,k]=capon2(x1,x2,M,METHOD);

% Default to foward-backward covariance estimation
if nargin < 4
    METHOD = 'modified';
end

N = length(x1);  % Signal length
L = N - M + 1;

R11 = corr_matlab(x1,x1,M,METHOD);
R22 = corr_matlab(x2,x2,M,METHOD);
R12 = corr_matlab(x1,x2,M,METHOD);
Ri11 = inv(R11);
Ri22 = inv(R22);

m = [0:M-1]';           %
k = [0:N/2]/N;        % Estimate positive and negative frequencies

% Make some pre-calculations to speed up things
Ri11_R12_Ri22 = Ri11*R12*Ri22;
Ri11_Ri22 = Ri11*Ri22;

% Calculate the cross PSD at each frequency
for l = 1:length(k)
   aM = exp(j*2*pi*k(l)*m);
   C(l) = (aM' * Ri11_R12_Ri22 * aM) / (aM' * Ri11_Ri22 * aM);
end

C = [fliplr(C(2:end))  C(1:end-1)]';
k = [-fliplr(k(2:end)) k(1:end-1)]';
