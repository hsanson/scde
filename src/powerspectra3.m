%
%   File:      powerspectra3.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17

%   Description:
%        Power Spectrum Density estimation using filter bank based methods.
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%        - The lse and APES methods do not work yet.

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 256;                        % Number of observed samples
A = [ 1 1 ];                    % Amplitudes of each sinusoid
f = [0.3;0.2];                  % Frequencies of each sinusoid
fs        = 1;                  % Sampling frequency
ts        = 1/fs;
snr       = 20;                 % Signal to noise ratio in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nT        = [0:N-1]*ts;             % Time axis

y = A*sin(2*pi*f*nT);               % generate sin signal
y = awgn(y,snr,'measured');         % Add awgn noise to the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation via LSE method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,kl] = lse(y,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation via WLS method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[W,kw] = wls(y,y,length(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation via CAPON method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C,kc] = capon(y,y,25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation via CAPON method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C2,k2] = capon2(y,y,8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation via APES method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ca,ka] = apes(y,y,length(nT)/4);

plot(kw*fs,10*log10(W),kc*fs,10*log10(C),k2*fs,10*log10(C2),ka*fs,10*log10(Ca));
legend('WLS','Capon','Capon2','APES');
title('Spectral Density');
ylabel('PSD db');
xlabel('Hz');

