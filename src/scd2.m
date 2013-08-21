%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      scd2.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/13

%   Description:
%       Estimates the Spectral Correlation Density of a signal using the
%       signal modulates method. This is similar to scd.m but here we calculate
%       the the SCD at a single cyclic frequency.
%
%   Notes:
%        - Tested with Matlab 2007a

function [SCD F] = scd(x,M,METHOD,fs,a)

% x      -  Process for which the SCD is to be estimated
% M      -  For WLS method means the segment size
%           For Capon/APES means the filter length
% METHOD -  Cross spectrum method = [ WLS, MLM, NMLM, APES ] 
% fs     -  Sampling frequency
% a      -  Cyclic frequency at wich to estimate SCD
%
% SCD    - The SCD at cyclic frequency a
% F      - The frequency axis in Hz

if exist('METHOD','var') == 0
    METHOD = 'wosa';
else
    %% Try to match one of the available methods
    methods_list = {'wosa','mvdr','nmvdr','apes'};
    indx = strmatch(lower(METHOD),methods_list);

    if isempty(indx) | length(indx) > 1,
        error('Ambiguous or invalid method specified.');
    end
    METHOD = methods_list{indx};
end

if exist('fs','var') == 0
    fs = 1;
end

x = x(:);  % make sure the vector is column ordered

ts = 1/fs;

N = length(x);      % Length of the process
n = [0:N-1]'*ts;

if strcmp(METHOD,'wosa')
    u = x .* exp(-j*pi*a*n);
    v = x .* exp(j*pi*a*n);
    [SCD f] = wls(u,v,M);
elseif strcmp(METHOD,'mvdr')
    u = x .* exp(-j*pi*a*n);
    v = x .* exp(j*pi*a*n);
    [SCD f] = capon(u,v,M);
elseif strcmp(METHOD,'nmvdr')
    u = x .* exp(-j*pi*a*n);
    v = x .* exp(j*pi*a*n);
    [SCD f] = capon2(u,v,M);
elseif strcmp(METHOD,'apes')
    u = x .* exp(-j*pi*a*n);
    v = x .* exp(j*pi*a*n);
    [SCD f] = apes(u,v,M);
else
    error(['Unknown method ' METHOD]); 
end


F = f*fs;   % Convert to Hz

