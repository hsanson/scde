%
%   File:      powerspectra7.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/24

%   Description:
%       Example of SCD estimation using different cross-spectrum estimators.
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14)

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = 16;            % Number of samples
fs = 1;             % Sampling Frequency
nsamp = 16;         % Pulse shaping oversampling
SNR = 20;           % SNR in DB

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
noise_var = sqrt(No/2)

noise = noise_var * randn(1, length(sn));

sn_noisy = sn + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theorical SCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [-Nup/2:Nup/2-1]*fs/Nup;
f = [-Nup/2:Nup/2-1]*fs/Nup;

S = zeros(length(f),length(a));
for i=1:length(a)

    Q1 = sin(pi*td*(f + fc + a(i)/2))./(pi*(f + fc + a(i)/2));
    Q2 = sin(pi*td*(f + fc - a(i)/2))./(pi*(f + fc - a(i)/2));
    Q3 = sin(pi*td*(f - fc + a(i)/2))./(pi*(f - fc + a(i)/2));
    Q4 = sin(pi*td*(f - fc - a(i)/2))./(pi*(f - fc - a(i)/2));

    S1(:,i) =  (1/(2*td)) * ( (Q1.*Q2 + Q3.*Q4) + (Q1.*Q4) + (Q2.*Q3) );
end

figure(3);
subplot(221);
contour(a,f,abs(S1));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('BPSK Theorical Cyclic Spectrum','FontSize', 16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the SCDs via different methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Sw fw aw] = scd(sn_noisy,N,'WLS',fs);
[Cs1 k1 a1] = scd(sn_noisy,N,'MLM',fs); 
[Cs2 k2 a2] = scd(sn_noisy,N,'NMLM',fs); 

subplot(222);
contour(aw,fw,abs(Sw));
title('BPSK Cyclic Spectrum (Welch)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(223);
contour(a1,k1,abs(Cs1));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('BPSK Cyclic Spectrum (CAPON-1)','FontSize', 16);

subplot(224);
contour(a2,k2,abs(Cs2));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('BPSK Cyclic Spectrum (CAPON-2)','FontSize', 16);

figure(4);

subplot(221);
mesh(a,f,abs(S1));
title('BPSK Theorical Cyclic Spectrum','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(222);
mesh(aw,fw,abs(Sw));
title('BPSK Cyclic Spectrum (Welch)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(223);
mesh(a1,k1,abs(Cs1));
title('BPSK Cyclic Spectrum (CAPON-1)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(224);
mesh(a2,k2,abs(Cs2));
title('BPSK Cyclic Spectrum (CAPON-2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);




% QPSK

%  res = zeros(length(f),length(a));
%
%  for i=1:length(a)
%
%      Q1 = sin(pi*td*(f + fc + a(i)/2))./(pi*(f + fc + a(i)/2));
%      Q2 = sin(pi*td*(f + fc - a(i)/2))./(pi*(f + fc - a(i)/2));
%      Q3 = sin(pi*td*(f - fc + a(i)/2))./(pi*(f - fc + a(i)/2));
%      Q4 = sin(pi*td*(f - fc - a(i)/2))./(pi*(f - fc - a(i)/2));
%
%      res(:,i) =  (1/(2*td)) * (Q1.*Q2 + Q3.*Q4);
%  end
%
%  figure(2);
%  mesh(res);
%  title('QPSK Theorical Cyclic Spectrum','FontSize', 24);


