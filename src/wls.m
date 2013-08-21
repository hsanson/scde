%
%   File:      wls.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/22

%   Description:
%       Obtains the cross Power Spectrum Density Using WLS method  of the processes
%       x1 and x2. For M = 1 we get the Periodogram and for M > 1 we get the averaged Welch
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%        - Based on the paper 
%           "Amplitude Estimation of Sinusoidal Signals: Survey, New Results, and an Application"
%           by Petre Stoica, Hongbin Li and Jian Li
%
%%

function [S,k]=wls(x1,x2,M,METHOD);

%  S,  cross PSD
%  x1, input sequence of length n
%  x2, input sequence of length n
%  M,  Filter length
%  METHOD, covariance estimation method

% Default to foward-backward covariance estimation
if nargin < 4
    METHOD = 'modified';   
end

N = length(x1);  % Signal length
L = N - M + 1;   %

R12 = corr_matlab(x1,x2,M,METHOD);

s = zeros(M,N);
k = [-N/2:N/2-1]/N;
for idx = 1:M
    s(idx,:) = exp(j*2*pi*(idx-1)*k);
end

S = (1/M) * real(diag(s'*R12*s));

%S = [fliplr(S(2:end)) ; S(1:end-1)]';  % Obtain the twosided version
%k = [-fliplr(k(2:end)) k(1:end-1)]'; 
 

%S = 10*log10((1/M) * real(diag(s'*R12*s)));  % dB scale
