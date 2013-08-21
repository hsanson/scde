%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      hankel_backward.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%     Builds the backward hankel matrix of vector x

function X=hankel_backward(x,M);

x = x(:)';

n = length(x);
l1 = 1;
l2 = M;
X = [];

while l2<=n
    X=[X;fliplr(x(l1:l2))];
    l1=l1+1;
    l2=l2+1;
end


