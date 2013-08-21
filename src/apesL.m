%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      apesL.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17
%
%   Description:
%       Estimates the cross PSD of two signals using the APES estimator.
%       This implementation is based on Dr. Lagunas code (CTTC) that creates a
%       fourier matrix and calculate the spectrum for all frequencies in one
%       step.
%
%   Notes:
%        - WARNING: This script does not work!! use apes.m instead.

function [S,k]=apes(x1,x2,M);

N = length(x1);  % Signal length
L = N - M + 1;   %

X1 = hankel_foward(x1,M);
X2 = hankel_foward(x2,M);

% Obtain the covariance matrices
R11 = corr_backward(x1,x1,M);
R22 = corr_backward(x2,x2,M);
R12 = corr_backward(x1,x2,M);

k = [-N:N-1]/(N);        % Estimate positive and negative frequencies

s = zeros(M,2*N);
k = [-N:N-1]/(2*N);
for idx = 1:M
    s(idx,:) = exp(j*2*pi*(idx-1)*k);
end

sl = zeros(L,2*N);
for idx = 1:L
    sl(idx,:) = exp(j*2*pi*(idx-1)*k);
end


gW1 = (1/L) * X1' * sl;
gW2 = (1/L) * X2' * sl;
Q11 = inv(R11 - gW1*gW1');
Q22 = inv(R22 - gW2*gW2');

num = s'*Q11*R12*Q22*s;
den = s'*Q11*s*s'*Q22*s;

S = (real(diag(num./den)));

