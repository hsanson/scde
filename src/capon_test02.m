%
%   File:      capon_test02.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/6
%
%   Description:
%       Example of SCD estimation of Sinusoid signal using MVDR for diferent filter sizes.
%   Notes:
%        - Tested with Matlab 2007a

clear; clc;
warning off all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = 64;              % Number of observed samples
A       = 1;               % SNR of each sinusoid
fc      = 1.5;             % Frequencies of each sinusoid
fs      = 10;              % Sampling frequency
trials  = 10;              % Number of monte carlo trials
snr     = 10;              % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
M       = [N/16 N/8 N/4 N/2]; % Filter lengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*sin(2*pi*fc*nT);
y = awgn(y,snr,'measured');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the SCDs via different methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,f,a] = scd(y,N/16,'MVDR',fs);
figure(1);
subplot(221);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(221);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/8,'MVDR',fs);
figure(1);
subplot(222);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(222);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (MVDR N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/4,'MVDR',fs);
figure(1);
subplot(223);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(223);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/2,'MVDR',fs);
figure(1);
subplot(224);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(224);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (MVDR N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

