%
%   File:      variance4.m
%   Author(s): Horacio Sanson

%   Description:
%     APES vs CAPON
%

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters

N  = 128;     % Number of samples
M  = [4 8 16 32 64 128 256 512 1024];
A = [1];                  % Amplitudes of each sinusoid
f = [0.07];      % Frequencies of each sinusoid
fs        = 5;                      % Sampling frequency
ts        = 1/fs;
variance  = 0.1;                    % Noise variance of the signal

trials = 10;                        % Number of trials

%% Calculate the spectra of the first value of N here so we can graph it

% Build the M vector
m = M(M<=N/2);

var_capon1 = zeros(size(m));
var_capon2 = zeros(size(m));
var_apes = zeros(size(m));
mean_capon1 = zeros(size(m));
mean_capon2 = zeros(size(m));
mean_apes = zeros(size(m));
rms_capon1 = zeros(size(m));
rms_capon2 = zeros(size(m));
rms_apes = zeros(size(m));
pwd_capon1 = zeros(size(m));
pwd_capon2 = zeros(size(m));
pwd_apes = zeros(size(m));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal generation

% Cool trick that generates a signal composed of the sum of
% several sinusoids + white noise. (Matlab Specral Analysis Help)
nT        = [0:N-1]*ts;              % Time axis
sn_noisy = A*cos(2*pi*f*nT) + variance * randn(size(nT));
pwd_theori = zeros(size(m)) + sum(abs(sn_noisy).^2)/length(sn_noisy);

Nfft = pow2(nextpow2(length(sn_noisy)));
nF = [-Nfft/2:Nfft/2-1]*fs/Nfft;          % Freq axis
scaling =  Nfft/fs;
    
for i = 1:length(m)

    %% Result vectors
    sum_capon1 = zeros(N,1);
    sum_capon2 = zeros(N,1);
    sum_apes = zeros(N,1);
    sqr_capon1 = zeros(N,1);
    sqr_capon2 = zeros(N,1);
    sqr_apes = zeros(N,1);
    
    for trial = 1:trials
        
        [Scapon1,k1] = capon(sn_noisy,sn_noisy,m(i),'modified');
        sum_capon1 = sum_capon1 + Scapon1;
        sqr_capon1 = sqr_capon1 + Scapon1.^2;
        pwd_capon1(i) = pwd_capon1(i) + sum(abs(real(Scapon1)).^2);
        [Scapon2,k2] = capon2(sn_noisy,sn_noisy,m(i),'modified');
        sum_capon2 = sum_capon2 + Scapon2;
        sqr_capon2 = sqr_capon2 + Scapon2.^2;
        pwd_capon2(i) = pwd_capon2(i) + sum(abs(real(Scapon2)).^2);
        [Sapes,k2] = apes(sn_noisy,sn_noisy,m(i),'modified');
        sum_apes = sum_apes + Sapes;
        sqr_apes = sqr_apes + Sapes.^2;
        pwd_apes(i) = pwd_apes(i) + sum(abs(real(Sapes)).^2);

    end

    sum_capon1 = sum_capon1./trials;                      %% Calculate the means per frequency 
    sqr_capon1 = sqr_capon1./trials - sum_capon1.*sum_capon1;    %% Calculate the variance per frequency 
    var_capon1(i) = mean(sqr_capon1) + var(sum_capon1);    %% Law of total variance or variance decomposition
    %var_capon1(i) = mean(sqr_capon1);    %% Law of total variance or variance decomposition
    mean_capon1(i) = mean(sum_capon1);
    rms_capon1(i) = sqrt(var_capon1(i).^2 + mean_capon1(i).^2);
    pwd_capon1(i) = pwd_capon1(i)./trials;
    
    sum_capon2 = sum_capon2./trials;                      %% Calculate the means per frequency 
    sqr_capon2 = sqr_capon2./trials - sum_capon2.*sum_capon2;    %% Calculate the variance per frequency 
    var_capon2(i) = mean(sqr_capon2) + var(sum_capon2);    %% Law of total variance or variance decomposition
    %var_capon2(i) = mean(sqr_capon2);
    mean_capon2(i) = mean(sum_capon2);
    rms_capon2(i) = sqrt(var_capon2(i).^2 + mean_capon2(i).^2);
    pwd_capon2(i) = pwd_capon2(i)./trials;

    sum_apes = sum_apes./trials;                      %% Calculate the means per frequency 
    sqr_apes = sqr_apes./trials - sum_apes.*sum_apes;    %% Calculate the variance per frequency 
    var_apes(i) = mean(sqr_apes) + var(sum_apes);    %% Law of total variance or variance decomposition
    %var_apes(i) = mean(sqr_apes);
    mean_apes(i) = mean(sum_apes);
    rms_apes(i) = sqrt(var_apes(i).^2 + mean_apes(i).^2);
    pwd_apes(i) = pwd_apes(i)./trials;
end

figure(2);
plot(m,var_capon1,'g*-');
hold on;
plot(m,var_capon2,'c+-');
plot(m,var_apes,'b+-');
title('PSD Variance','FontSize', 10);
xlabel('M (Filter Size)','FontSize', 8);
ylabel('Variance','FontSize', 8);
h = legend('Capon1','Capon2','Apes');
set(h,'FontSize',10,'Location','NorthEast');

figure(3);
plot(m,mean_capon1,'g*-');
hold on;
plot(m,mean_capon2,'c+-');
plot(m,mean_apes,'b+-');
title('PSD Mean','FontSize', 10);
xlabel('M (Filter Size)','FontSize', 8);
ylabel('Mean','FontSize', 8);
h = legend('Capon1','Capon2','Apes');
set(h,'FontSize',10,'Location','NorthEast');

figure(4);
plot(m,rms_capon1,'g*-');
hold on;
plot(m,rms_capon2,'c+-');
plot(m,rms_apes,'b+-');
title('PSD RMS','FontSize', 10);
xlabel('M (Filter Size)','FontSize', 8);
ylabel('RMS','FontSize', 8);
h = legend('Capon1','Capon2','Apes');
set(h,'FontSize',10,'Location','NorthEast');

figure(5);
plot(m,pwd_capon1,'g*-');
hold on;
plot(m,pwd_capon2,'c+-');
plot(m,pwd_apes,'b+-');
plot(m,pwd_theori,'r:');
title('Average Power','FontSize', 10);
xlabel('M (Filter Size)','FontSize', 8);
ylabel('RMS','FontSize', 8);
h = legend('Capon1','Capon2','Apes','Theorical');
set(h,'FontSize',10,'Location','NorthEast');


