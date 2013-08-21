%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      lse.m
%   Author(s): Horacio Sanson
%
%   Description:
%       Estimates the cross PSD of two signals using the LSE estimator.
%
function [S,k]=lse(x1,x2,METHOD);

% Default to foward-backward covariance estimation
if nargin < 4
    METHOD = 'modified';   
end

N = length(x1);  % Signal length

R12 = corr_matlab(x1,x2,N,METHOD);
k = [-N/2:N/2-1]/N;
s = exp(-j*2*pi*k');

num = s'*R12*s;
den = s'*s;

%S = 10*log10(real(num/den));
S = num/den;

