%
%   File:      variance1.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/18

%   Description:
%     Periodogram vs Welch RMS
%     Note that RMS^2 = mean^2 + var^2

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters

N  = [32:20:512]';      % Number of samples
A = [ 1 1 ];            % Amplitudes of each sinusoid
f = [0.08;0.2];         % Frequencies of each sinusoid
fs        = 1;          % Sampling frequency
ts        = 1/fs;
variance  = 0.1;        % Noise variance of the signal

trials = 10;            % Number of trials

var_welch = zeros(size(N));
var_capon1 = zeros(size(N));
var_capon2 = zeros(size(N));
mean_welch = zeros(size(N));
mean_capon1 = zeros(size(N));
mean_capon2 = zeros(size(N));

for i = 1:length(N)

    nT        = [0:N(i)-1]*ts;      % Time axis
    Sperio = zeros(trials,N(i));
    Swelch = zeros(trials,N(i));
    Scapon1 = zeros(trials,N(i));
    Scapon2 = zeros(trials,N(i));
    Nfft = N(i);
    M = min(round(N(i)/4),25);      % Filter length. Make sure it does not exceed 25

    for trial=1:trials

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Signal generation

        % Cool trick that generates a signal composed of the sum of
        % several sinusoids + white noise. (Matlab Specral Analysis Help)
        sn_noisy = A*cos(2*pi*f*nT) + variance * randn(size(nT));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [Swelch(trial,:),nF] =  wls(sn_noisy,sn_noisy,N(i));

        [Scapon1(trial,:),k1] = capon(sn_noisy,sn_noisy,M);

        [Scapon2(trial,:),k2] = capon2(sn_noisy,sn_noisy,M);


    end

    var_welch(i) = sum(var(Swelch,1));
    var_capon1(i) = sum(var(Scapon1,1));
    var_capon2(i) = sum(var(Scapon2,1));
    mean_welch(i) = sum(var(Swelch,1));
    mean_capon1(i) = sum(var(Scapon1,1));
    mean_capon2(i) = sum(var(Scapon2,1));

    %figure(i);
    %plot(nF*fs,10*log10(Swelch'),'r-',k1*fs,10*log10(real(Scapon1')),'g-',k2*fs,10*log10(real(Scapon2')),'c-');
    %title(['Estimated PSD (N=' int2str(N(i)) ')'],'FontSize', 10);
    %xlabel('Hz','FontSize', 8);
    %ylabel('PSD dB','FontSize', 8);
    %h = legend('WLS','Capon1','Capon2');
    %set(h,'FontSize',10,'Location','SouthEast');
    %hold off;
end


figure(10);
plot(N,var_welch,'r-');
hold on;
plot(N,var_capon1,'g*-');
plot(N,var_capon2,'c+-');
title('PSD Variance','FontSize', 10);
xlabel('N (Number of samples)','FontSize', 8);
ylabel('Variance','FontSize', 8);
h = legend('WLS','Capon1','Capon2');
set(h,'FontSize',10,'Location','NorthEast');

figure(11);
plot(N,mean_welch,'r-');
hold on;
plot(N,mean_capon1,'g*-');
plot(N,mean_capon2,'c+-');
title('PSD Mean','FontSize', 10);
xlabel('N (Number of samples)','FontSize', 8);
ylabel('Mean','FontSize', 8);
h = legend('WLS','Capon1','Capon2');
set(h,'FontSize',10,'Location','NorthEast');

%figure(12);
%plot(N,rms_welch,'r-');
%hold on;
%plot(N,rms_capon1,'g*-');
%plot(N,rms_capon2,'c+-');
%title('PSD RMS','FontSize', 10);
%xlabel('N (Number of samples)','FontSize', 8);
%ylabel('RMS','FontSize', 8);
%h = legend('Welch','Capon1','Capon2');
%set(h,'FontSize',10,'Location','NorthEast');



