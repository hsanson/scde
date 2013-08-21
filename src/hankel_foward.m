%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      hankel_backward.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%     Builds the foward M x L hankel matrix of vector x
%     where L = length(x) - M + 1
%
%          x(0)    x(1)   . . .  x(L-2)   x(L-1)
%     X =  x(1)    x(2)   . . .  x(L-1)   x(L)
%          :       :             :        :
%          x(M-2)  x(M-1) . . .  x(N-3)   x(N-2)
%          x(M-1)  x(M)   . . .  x(N-2)   x(N-1)
%

function X=hankel_foward(x,M);

x = x(:)';

n = length(x);
l1 = 1;
l2 = M;
X = [];

while l2<=n
    X=[X;x(l1:l2)];
    l1=l1+1;
    l2=l2+1;
end

X = X';

