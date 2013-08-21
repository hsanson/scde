%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      bpsk_test1.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%     BPSK PSD estimation test
clear; clc;

M = 2;                              % Alphabet Size
k = log2(M);                        % Number of Bits per symbol
fs = 1;                             % Sampling frequency
T = 1/fs;                           % Sampling period
t = [0:T:100];                      % time vector
nsamp = 2;                          % Oversampling rate
fftlength = 2^nextpow2(length(t));  % FFT length
ebno = 25;                           % Eb/No

% Initialize the root-raised cosine (Pulse Shaping Filter)
filtorder = 20;                 % Filter order
delay = filtorder/(nsamp*2);    % Group delay (# of input samples)
rolloff = 0.3;                  % Rolloff factor of filter

% Create a square root raised cosine filter.
rrcfilter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

% Signal Source
x = randint(length(t),1);   % Random binary data

% Now transmit the signal throught a channel with
% different EbNo values.

% Modulate the signal
ytx = pskmod(x,M);

% Pulse Shape and oversample the signal
ytx = rcosflt(ytx,1,nsamp,'filter',rrcfilter);

% add some AWGN noise
snr = ebno + 10*log10(k) - 10*log10(nsamp);
ynoisy = awgn(ytx, snr, 'measured');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theorical BPSK PSD (Raised Cosine)
%% from Proakis - Digital Communications 4th edition pg 547
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [0:fftlength/2]*fs/fftlength;
Px = T/2 * (1 + cos((pi*T/rolloff)*(abs(f)-((1 - rolloff)/(2*T)))));
Px(f < (1-rolloff)/(2*T)) = T;
Px(f > (1+rolloff)/(2*T)) = 0;

Px = var(ynoisy)*(abs(Px).^2)./T;

figure(1);
grid on;
subplot(311);
plot(f,10*log10(Px/max(Px)))
title('Theorical BPSK PSD');
xlabel('Hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSD estimation via Matlab periodogram method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Px = pwelch(ynoisy, [],[], fftlength, fs);
subplot(312);
f = [0:fftlength/2]*fs/fftlength;
plot(f,10*log10(Px(1:length(f))/max(Px)), 'r');
title('Periodogram Method');
xlabel('Hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSD estimation via Periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = fft(ynoisy,fftlength);
f = [0:fftlength/2]*fs/fftlength;
Px=X.*conj(X)./length(ynoisy);

subplot(313);
plot(f,10*log10(Px(1:length(f))), 'r');
title('Estimated BPSK PSD (FFT method)');
xlabel('Hz');


