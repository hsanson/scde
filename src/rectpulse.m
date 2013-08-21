%
%   File:      rectpulse.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/18

%   Description:
%       RECTPULSE Rectangular pulse shaping (Mimick Matlab function of same name).
%       Y = RECTPULSE(X,N) returns Y, a rectangular pulse shaped version of X,
%       with N samples per symbol. This function replicates each symbol in
%       X N times.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9a
%        - Based on the octave-forge implementation of upsample.m by Paul Kienzle

function y = rectpulse(x,n)

%Check x, n
if( ~isnumeric(x))
    error('rectpulse:x','X must be numeric.');
end

if(~isreal(n) || ~isscalar(n) ||  n<=0 || (ceil(n)~=n) || ~isnumeric(n) )
    error('rectpulse:n','N must be a positive integer.');
end

[nr,nc] = size(x);
if any([nc,nr]==1)
    y = zeros(n*nr*nc,1);
    for phase = 1:n
        y(phase:n:end) = x;
    end
    if nr==1, y=y.'; end
else
    y = zeros(n*nr,nc);
    for phase = 1:n
        y(phase:n:end,:) = x;
    end
end

