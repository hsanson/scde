%
%   File:      corr_matlab.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17

%   Description:
%       Obtains the MxM cross covariance matrix of vectors x and y with M levels of freedom.
%       This is a simple wrapper on top of the corrmtx function to get the cross covariance matrix
%       of two signals.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%
%
function R=corr_matlab(x,y,M,METHOD);

if nargin < 4
    METHOD = 'modified';   % Use foward-backward method by default
end

if(M <= 0) 
    error('M <= 0 calculating correlation matrix');
end

% Obtain the covariance matrices
A1 = corrmtx(x,M-1,METHOD);
A2 = corrmtx(y,M-1,METHOD);
R = (A1'*A2 + (A1'*A2)')/2;

