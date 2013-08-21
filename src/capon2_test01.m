%
%   File:      capon2_test01.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/6
%
%   Description:
%       Example of SCD estimation using NMVDR for diferent filter sizes
%   Notes:
%        - Tested with Matlab 2007a

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = 16;        % Number of samples
fs = 1;         % Sampling Frequency
nsamp = 16;     % Pulse shaping oversampling
SNR = 20;       % SNR in DB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts = 1/fs;
fd = 1/(8*ts);
td = 1/fd;
fc = 1/(5*ts);
snr = 10^(SNR/10);          % SNR as ratio

Nup = N*nsamp;      % Number of samples of the upsampled signal

t = [0:Nup-1]*ts;     % time vector
Nfft = pow2(nextpow2(Nup));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPSK Signal generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn = sign(rand(1,N)-0.5);   % Create BPSK signal
sn = rectpulse(sn,nsamp);

win = rectwin(Nfft)';

%% Carrier modulate
sn = sn.*cos(2*pi*fc*t);

%% Add some Gaussian noise
Eb = sum(abs(sn).^2)/length(sn);
No = Eb/snr;
noise_var = sqrt(No/2);

noise = noise_var * randn(1, length(sn));

sn_noisy = sn + noise;

L = length(sn_noisy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the SCDs via different methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,f,a] = scd(sn_noisy,L/16,'NMVDR',fs);
figure(1);
subplot(221);
contour(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(221);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(sn_noisy,L/8,'NMVDR',fs);
figure(1);
subplot(222);
contour(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(222);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(sn_noisy,L/4,'NMVDR',fs);
figure(1);
subplot(223);
contour(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(223);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(sn_noisy,L/2,'NMVDR',fs);
figure(1);
subplot(224);
contour(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(224);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (NMVDR N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

