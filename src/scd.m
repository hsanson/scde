%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      scd.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/1

%   Description:
%       Estimates the Spectral Correlation Density of a signal using the
%       signal modulates method.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%        - Tested with Matlab 2007a

function [SCD F A] = scd(x,M,METHOD,fs,res)

% x      -  Process for which the SCD is to be estimated
% M      -  For WLS method means the segment size
%           For Capon/APES means the filter length
% METHOD -  Cross spectrum method = [ WLS, MLM, NMLM, APES ] 
% fs     -  Sampling frequency
% res    -  cycle resolution
%
% SCD    - The bi-frequency SCD
% F      - The frequency axis in Hz
% A      - The cyclic frequency axis in Hz

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

if exist('res','var') == 0
    res = length(x);
end

x = x(:);  % make sure the vector is column ordered

ts = 1/fs;

N = length(x);      % Length of the process
n = [0:N-1]'*ts;
A = [0:res/2]*fs/res;   % Obtain the cyclic frequency axis from 0 to N or 0 to fs

if strcmp(METHOD,'wosa')
    for i = 1:length(A)
        u = x .* exp(-j*pi*A(i)*n);
        v = x .* exp(j*pi*A(i)*n);
        [SCD(:,i) f] = wls(u,v,M);
    end
elseif strcmp(METHOD,'mvdr')
    for i = 1:length(A)
        u = x .* exp(-j*pi*A(i)*n);
        v = x .* exp(j*pi*A(i)*n);
        [SCD(:,i) f] = capon(u,v,M);
    end
elseif strcmp(METHOD,'nmvdr')
    for i = 1:length(A)
        u = x .* exp(-j*pi*A(i)*n);
        v = x .* exp(j*pi*A(i)*n);
        [SCD(:,i) f] = capon2(u,v,M);
    end
elseif strcmp(METHOD,'apes')
    for i = 1:length(A)
        u = x .* exp(-j*pi*A(i)*n);
        v = x .* exp(j*pi*A(i)*n);
        [SCD(:,i) f] = apes(u,v,M);
    end
else
    error(['Unknown method ' METHOD]); 
end


F = f*fs;   % Convert to Hz

%% We calculate the SCD for cyclic frequencies from 0 to fs/2
%% so now following the symetry property we obtain the SCD for
%% negative frequencies by mirroring.

SCD = [fliplr(conj(SCD(:,2:end))) SCD(:,1:end-1)];
A = [-fliplr(A(:,2:end)) A(:,1:end-1)];


