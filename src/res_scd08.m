% 
%   File:      res_scd04.m
%   Author(s): Horacio Sanson      
%   Revision : 2007/11/14

%   Description: 
%        Comparison of probability of resolution of SCD estimators for
%        different SNR and filter lengths and plots them as functions 
%        of varying SNR
%
%   Notes: 
%       - Tested with Matlab 2007a
%       - Resolution criterion S(wm) = 1/2(S(w1) + S(w2)) taken from the paper
%         Probability of Resolution of the MUSIC algorithm


clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = 64;              % Number of observed samples
A       = [1 1];           % Amplitude of each sinusoid
fc      = [1.0 1.5]';      % Frequencies of each sinusoid
fs      = 10;              % Sampling frequency
trials  = 1;             % Number of monte carlo trials
snr     = [-100:1:100];      % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
a       = 2*fc;               % cyclic frequency two

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*cos(2*pi*fc*nT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pres = zeros(length(snr),5);   % probability of resolution

for idx2 = 1:length(snr)

    success_wosa = 0;
    success_mvdr = 0;
    success_nmvdr = 0;
    success_nmvdr2 = 0;
    success_apes = 0;

    for trial=[1:trials]    % Repeat experiment trials times

        disp([' SNR ' int2str(snr(idx2)) ' trial ' int2str(trial)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add gausian noise to the signals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        yn = awgn(y,snr(idx2),'measured');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SCD Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [S1 f1] = scd2(yn,N/4,'WOSA',fs,2*fc(1));
        [S2 f2] = scd2(yn,N/4,'WOSA',fs,2*fc(2));
        [S3 f3] = scd2(yn,N/4,'WOSA',fs,fc(1)+fc(2));

        % Find the bin at frequency 0
        [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
        [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
        [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

        if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
            success_wosa = success_wosa + 1;
        end

        [S1 f1] = scd2(yn,N/4,'MVDR',fs,2*fc(1));
        [S2 f2] = scd2(yn,N/4,'MVDR',fs,2*fc(2));
        [S3 f3] = scd2(yn,N/4,'MVDR',fs,fc(1)+fc(2));

        % Find the bin at frequency 0
        [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
        [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
        [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

        if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
            success_mvdr = success_mvdr + 1;
        end

        [S1 f1] = scd2(yn,N/4,'NMVDR',fs,2*fc(1));
        [S2 f2] = scd2(yn,N/4,'NMVDR',fs,2*fc(2));
        [S3 f3] = scd2(yn,N/4,'NMVDR',fs,fc(1)+fc(2));

        % Find the bin at frequency 0
        [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
        [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
        [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

        if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
            success_nmvdr = success_nmvdr + 1;
        end

        [S1 f1] = scd2(yn,N/8,'NMVDR',fs,2*fc(1));
        [S2 f2] = scd2(yn,N/8,'NMVDR',fs,2*fc(2));
        [S3 f3] = scd2(yn,N/8,'NMVDR',fs,fc(1)+fc(2));

        % Find the bin at frequency 0
        [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
        [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
        [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

        if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
            success_nmvdr2 = success_nmvdr2 + 1;
        end

        [S1 f1] = scd2(yn,N/4,'APES',fs,2*fc(1));
        [S2 f2] = scd2(yn,N/4,'APES',fs,2*fc(2));
        [S3 f3] = scd2(yn,N/4,'APES',fs,fc(1)+fc(2));

        % Find the bin at frequency 0
        [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
        [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
        [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

        if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
            success_apes = success_apes + 1;
        end

    end

    pres(idx2,1)  = success_wosa/trials;
    pres(idx2,2)  = success_mvdr/trials;
    pres(idx2,3)  = success_nmvdr/trials;
    pres(idx2,4)  = success_nmvdr2/trials;
    pres(idx2,5)  = success_apes/trials;
end

save 'res_scd08.mat';
%pres = load 'res_scd08.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation probability of resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,pres(:,1),'o-',snr,pres(:,2),'^-',snr,pres(:,3),'+-',snr,pres(:,4),'x-',snr,pres(:,5),'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N/4 length', 'MVDR N/4 length', 'NMVDR N/4 length', 'NMVDR N/8 length', 'APES N/4 length');
title('Probability of resolution of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation probability of resolution without APES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(2);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,pres(:,1),'o-',snr,pres(:,2),'^-',snr,pres(:,3),'+-',snr,pres(:,4),'x-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N/4 length', 'MVDR N/4 length', 'NMVDR N/4 length','NMVDR N/8 length');
title('Probability of resolution of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res07
