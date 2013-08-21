%
%   File:      powerspectra4.m
%   Author(s): Horacio Sanson

%   Description:
%
%     PSD of a BPSK signal
%
%     The spectrum of a baseband BPSK signal, if the symbols are
%     equiproblables and uncorrelated, is the same as the spectrum of the
%     shaping pulse used. In this example we use a rectangular pulse shape 
%     so the resulting spectrum resembles a sinc function.
%
%     Depending on the Upsampling (NSAMP) and the number of samples (N) the
%     spectrum more closely approximates the sinc function.

clear; clc;

% Simulation Parameters
fs      = 0.5;                           % Sampling Frequency (> 2fc)
ts      = 1/fs;                          % Sampling Period
fc      = 0.1;                           % Carrier Frequency
M       = 2;                             % Alphabet Size (BPSK - 2. QPSK - 4)
N       = 16;                            % Number of samples
ebno    = 30;                            % EbNo in DB
NSAMP   = 8;                             % Upsampling?
trials  = 100;                           % Number of trials to run the estimation

% Derived parameters
k = log2(M);                            % Number of Bits per symbol

snr = ebno + 10*log10(k) - 10*log10(NSAMP); % SNR as ratio

% Do the estimation several times to see the variance
% of the estimates.
for i=1:trials
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal Source
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = randint(N,1);               % Random binary data

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BPSK Modulate the signal (M = 2) using Matlab functions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = pskmod(x,M);                        % Complex envelop??
    x = rectpulse(x,NSAMP);                 % Pulse Shaping
    L = length(x);
    n = [0:L-1]*ts;                 % Time Axis of pulse shaped signal
    %x = x.*cos(2*pi*fc.*n');               % Carrier modulation
    x = awgn(x, snr, 'measured');           %  Add AWGN noise


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PSD estimation via WLS method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Px f] = wls(x,x,L);
    plot(f*fs,10*log10(real(Px)), 'r');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation via CAPON method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on;
    [S,k] = capon(x,x,L/8);
    plot(k*fs,10*log10(real(S)), 'LineWidth', 2, 'Color', [0.6 0.7 1]);
end
title(['Estimated BPSK PSD over ' int2str(trials) ' iterations']);
xlabel('Hz');
ylabel('PSD (dB)');
legend('WLS method', 'Capon (MLM)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOISE http://www.dsprelated.com/showmessage/45398/1.php
