%
%   File:      variance3.m
%   Author(s): Horacio Sanson

%   Description:
%     Periodogram vs Welch RMS
%     Note that RMS^2 = mean^2 + var^2

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters

N  = [32 64 128 256 512 1024]';     % Number of samples
A = [ 1 1 1 1 1 ];                  % Amplitudes of each sinusoid
A = [ 1 ];                  % Amplitudes of each sinusoid
f = [0.05;0.06;0.07;0.08;0.2];      % Frequencies of each sinusoid
f = [0.07];      % Frequencies of each sinusoid
fs        = 5;                      % Sampling frequency
ts        = 1/fs;
variance  = 0.1;                    % Noise variance of the signal

trials = 100;                        % Number of trials



%% Calculate the spectra of the first value of N here so we can graph it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal generation

% Cool trick that generates a signal composed of the sum of
% several sinusoids + white noise. (Matlab Specral Analysis Help)
nT        = [0:N(1)-1]*ts;              % Time axis
sn_noisy = A*cos(2*pi*f*nT) + variance * randn(size(nT));
sumsqr(sn_noisy)/length(sn_noisy)

Nfft = pow2(nextpow2(length(sn_noisy)));
nF = [-Nfft/2:Nfft/2-1]*fs/Nfft;          % Freq axis
scaling =  Nfft/fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD estimation

%% PSD obtained via Periodogram
Sperio = fft(sn_noisy,Nfft);
Sperio = Sperio.*conj(Sperio)/Nfft;
sum(Sperio)
Sperio = 10*log10(Sperio/max(Sperio));
Sperio = fftshift(Sperio);

%% PSD obtained via Welch
Swelch =  pwelch(sn_noisy,[],[],Nfft,fs,'twosided');
Swelch = Swelch/scaling;
sum(Swelch)
Swelch = 10*log10(Swelch/max(Swelch));
Swelch = fftshift(Swelch);

%% PSD obtained via Capon-1 and Capon-2
Spmtm = pmtm(sn_noisy,[],Nfft,fs,'twosided');
Spmtm = Spmtm/scaling;
sum(Spmtm)
Spmtm = 10*log10(Spmtm/max(Spmtm));
Spmtm = fftshift(Spmtm);

figure(1)
plot(nF,Sperio,'b-');
hold on;
plot(nF,Swelch,'r+-');
plot(nF,Spmtm,'g+-');
title('Estimated PSD (N=32)','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);
h = legend('Periodogram','Welch','Multi Taper');
set(h,'FontSize',10,'Location','SouthEast');

a(0)
    
var_perio = zeros(size(N));
var_welch = zeros(size(N));
var_pmtm = zeros(size(N));
mean_perio = zeros(size(N));
mean_welch = zeros(size(N));
mean_pmtm = zeros(size(N));
rms_perio = zeros(size(N));
rms_welch = zeros(size(N));
rms_pmtm = zeros(size(N));

for i = 1:length(N)

    %% Result vectors
    sum_perio = zeros(N(i),1);
    sum_welch = zeros(N(i),1);
    sum_pmtm = zeros(N(i),1);
    sqr_perio = zeros(N(i),1);
    sqr_welch = zeros(N(i),1);
    sqr_pmtm = zeros(N(i),1);

    for trial=1:trials

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Signal generation

        % Cool trick that generates a signal composed of the sum of
        % several sinusoids + white noise. (Matlab Specral Analysis Help)
        nT        = [0:N(i)-1]*ts;              % Time axis
        sn_noisy = A*cos(2*pi*f*nT) + variance * randn(size(nT));

        Nfft = pow2(nextpow2(length(sn_noisy)));
        nF = [-Nfft/2:Nfft/2-1]*fs/Nfft;          % Freq axis
        scaling =  Nfft/fs;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PSD estimation

        %% PSD obtained via Periodogram
        Sperio =  periodogram(sn_noisy,[],Nfft,fs,'twosided');
        Sperio = Sperio/scaling;  % Matlab missing scale factor
        %Sperio = fft(sn_noisy.*hanning(Nfft)',Nfft)';
        %Sperio = (Sperio.*conj(Sperio))/Nfft;
        sum_perio = sum_perio + Sperio;
        sqr_perio = sqr_perio + Sperio.^2;
        Sperio = 10*log10(Sperio/max(Sperio));

        %% PSD obtained via Welch
        Swelch =  pwelch(sn_noisy,[],[],Nfft,fs,'twosided');
        Swelch = Swelch/scaling;  % Matlab missing scale factor
        sum_welch = sum_welch + Swelch;
        sqr_welch = sqr_welch + Swelch.^2;
        Swelch = 10*log10(Swelch/max(Swelch));

        Spmtm = pmtm(sn_noisy,[],Nfft,fs,'twosided');
        Spmtm = Spmtm/scaling;
        sum_pmtm = sum_pmtm + Spmtm;
        sqr_pmtm = sqr_pmtm + Spmtm.^2;
        Spmtm = 10*log10(Spmtm/max(Spmtm));

    end

    
    sum_perio = sum_perio./trials;                      %% Calculate the means per frequency 
    sqr_perio = sqr_perio./trials - sum_perio.*sum_perio;    %% Calculate the variance per frequency 
    var_perio(i) = mean(sqr_perio) + var(sum_perio);    %% Law of total variance or variance decomposition
    %var_perio(i) = sum(sqr_perio);
    mean_perio(i) = mean(sum_perio);
    rms_perio = sqrt(var_perio.^2 + mean_perio.^2);

    sum_welch = sum_welch./trials;                      %% Calculate the means per frequency 
    sqr_welch = sqr_welch./trials - sum_welch.*sum_welch;    %% Calculate the variance per frequency 
    var_welch(i) = mean(sqr_welch) + var(sum_welch);    %% Law of total variance or variance decomposition
    %var_welch(i) = sum(sqr_welch);
    mean_welch(i) = mean(sum_welch);
    rms_welch = sqrt(var_welch.^2 + mean_welch.^2);

    sum_pmtm = sum_pmtm./trials;                      %% Calculate the means per frequency 
    sqr_pmtm = sqr_pmtm./trials - sum_pmtm.*sum_pmtm;    %% Calculate the variance per frequency 
    var_pmtm(i) = mean(sqr_pmtm) + var(sum_pmtm);    %% Law of total variance or variance decomposition
    %var_pmtm(i) = sum(sqr_pmtm);    %% Law of total variance or variance decomposition
    mean_pmtm(i) = mean(sum_pmtm);
    rms_pmtm = sqrt(var_pmtm.^2 + mean_pmtm.^2);
    
end

figure(2);
plot(N,var_perio,'b-*');
hold on;
plot(N,var_welch,'r-');
plot(N,var_pmtm,'g*-');
title('PSD Variance','FontSize', 10);
xlabel('N (Number of samples)','FontSize', 8);
ylabel('Variance','FontSize', 8);
h = legend('Periodogram','Welch','Multi Taper');
set(h,'FontSize',10,'Location','NorthEast');

figure(3);
plot(N,mean_perio,'b-*');
hold on;
plot(N,mean_welch,'r-');
plot(N,mean_pmtm,'g*-');
title('PSD Mean','FontSize', 10);
xlabel('N (Number of samples)','FontSize', 8);
ylabel('Mean','FontSize', 8);
h = legend('Periodogram','Welch','Multi Taper');
set(h,'FontSize',10,'Location','NorthEast');

figure(4);
plot(N,rms_perio,'b-*');
hold on;
plot(N,rms_welch,'r-');
plot(N,rms_pmtm,'g*-');
title('PSD RMS','FontSize', 10);
xlabel('N (Number of samples)','FontSize', 8);
ylabel('RMS','FontSize', 8);
h = legend('Periodogram','Welch','Multi Taper');
set(h,'FontSize',10,'Location','NorthEast');

figure(5)
plot(nF,fftshift(Sperio),'b-');
hold on;
plot(nF,fftshift(Swelch),'r+-');
plot(nF,fftshift(Spmtm),'g+-');
title('Estimated PSD (N=1024)','FontSize', 10);
xlabel('Hz','FontSize', 8);
ylabel('PSD dB','FontSize', 8);
h = legend('Periodogram','Welch','Multi Taper');
set(h,'FontSize',10,'Location','SouthEast');



