%
%   File:      powerspectra5.m
%   Author(s): Horacio Sanson

%   Description:
%     Plot cross PSD of the BPSK signal using Capon-1
%     and Capon-2 with different values of filter length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = 16;            % Number of samples
fs = 1;           % Sampling Frequency
nsamp = 16;         % Pulse shaping oversampling
SNR = 20;           % SNR in DB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts = 1/fs;
fd = 1/(8*ts);
td = 1/fd;
fc = 1/(4*ts);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate and plot the PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
[C,k] = capon(sn_noisy,sn_noisy,length(sn_noisy)/16);
subplot(321);
plot(k*fs,C, 'b');
title('Capon-1 Method with M=N/16','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);

[C,k] = capon2(sn_noisy,sn_noisy,length(sn_noisy)/16);
subplot(322);
plot(k*fs,C, 'b');
title('Capon-2 Method with M=N/16','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);


[C,k] = capon(sn_noisy,sn_noisy,length(sn_noisy)/8);
subplot(323);
plot(k*fs,C, 'b');
title('Capon-1 Method with M=N/8','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);

[C,k] = capon2(sn_noisy,sn_noisy,length(sn_noisy)/8);
subplot(324);
plot(k*fs,C, 'b');
title('Capon-2 Method with M=N/8','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);

[C,k] = capon(sn_noisy,sn_noisy,length(sn_noisy)/4);
subplot(325);
plot(k*fs,C, 'b');
title('Capon-1 Method with M=N/4','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);

[C,k] = capon2(sn_noisy,sn_noisy,length(sn_noisy)/4);
subplot(326);
plot(k*fs,C, 'b');
title('Capon-2 Method with M=N/4','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);

