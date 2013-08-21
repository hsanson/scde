%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      apes.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17
%
%   Description:
%       Estimates the cross PSD of two signals using the APES estimator.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%        - We calculate only PSD for positive frequencies and mirror them to get the negative ones.
%        - WARNING: this script can take a LOT of time for filter sizes larger than 25.

function [C,k]=apes(x1,x2,M,METHOD);

% Default to foward-backward covariance estimation
if nargin < 4
    METHOD = 'modified';
end

N = length(x1);  % Signal length
L = N - M + 1;   %

X1 = hankel_foward(x1,M);
X2 = hankel_foward(x2,M);

% Obtain the covariance matrices
R11 = corr_matlab(x1,x1,M,METHOD);
R22 = corr_matlab(x2,x2,M,METHOD);
R12 = corr_matlab(x1,x2,M,METHOD);


l = [0:L-1]';
m = [0:M-1]';
k = [0:N/2]/N;                     % Estimate only positive frequencies

% Calculate the cross PSD
for ii = 1:length(k)
   aL = exp(-j*2*pi*k(ii)*l);
   aM = exp(j*2*pi*k(ii)*m);

   gW1 = (1/L) * X1 * aL;
   gW2 = (1/L) * X2 * aL;

   Q11 = inv(R11 - gW1*gW1');
   Q22 = inv(R22 - gW2*gW2');

   C(ii) = (aM' * Q11 * R12 * Q22 * aM) / (aM' * Q11 * aM * aM' * Q22 * aM);
end

C = [fliplr(C(2:end))  C(1:end-1)]';
k = [-fliplr(k(2:end)) k(1:end-1)]';

