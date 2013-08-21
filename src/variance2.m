%
%   File:      variance2.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/18

%   Description:
%       Obtain the variance of the PSD estimators for increasing number filter lengths
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%        - Here we show an example of how to obtain noise variance from the Eb/No value.
%          This would be the same as using the built in awgn function.

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 16;                     % Number of samples
nsamp = 16;                 % Pulse shaping oversampling
M = [4:4:32]';              % Filter lengths
fs = 1;                     % Sampling Frequency
SNR = 20;                   % SNR in DB
trials = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts = 1/fs;
fd = 1/(8*ts);
td = 1/fd;
fc = 1/(4*ts);
snr = 10^(SNR/10);          % SNR as ratio
len = N*nsamp;

var_capon1 = zeros(size(M));
var_capon2 = zeros(size(M));


for i = 1:length(M)

    psd_capon1 = zeros(trials,len);
    psd_capon2 = zeros(trials,len);

    for trial=1:trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BPSK Signal generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sn = sign(rand(1,N)-0.5);   % Create BPSK signal
        sn = rectpulse(sn,nsamp);

        t=[0:length(sn)-1];         % time vector
        Nfft = pow2(nextpow2(length(t)));

        %% Carrier modulate
        sn = sn.*cos(2*pi*fc*t);

        %% Add some Gaussian noise
        Eb = sum(abs(sn).^2)/length(sn);
        No = Eb/snr;
        noise_var = sqrt(No/2);

        noise = noise_var * randn(1, length(sn));

        sn_noisy = sn + noise;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PSD estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        psd_capon1(trial,:) = capon(sn_noisy,sn_noisy,M(i))';
        psd_capon2(trial,:) = capon2(sn_noisy,sn_noisy,M(i))';

        %Scapon1 = 10*log10(Scapon1/max(Scapon1));
        %Scapon2 = 10*log10(Scapon2/max(Scapon2));
    end

    var_capon1(i) = sum(var(psd_capon1));
    var_capon2(i) = sum(var(psd_capon2));

end

figure(1);
plot(M,var_capon1,'r+-',M,var_capon2,'g^-');
title('PSD Variance');
xlabel('length/M (Filter Length)');
ylabel('Variance');
legend('Capon1','Capon2');
