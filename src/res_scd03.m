% 
%   File:      res_scd03.m
%   Author(s): Horacio Sanson      
%   Revision : 2007/11/13

%   Description: 
%        Obtain the probability of resolution of NMVDR SCD estimators for different SNR
%        and filter lengths and plots them as functions of varying SNR
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
fc      = [1.0 2.0]';  % Frequencies of each sinusoid
fs      = 10;              % Sampling frequency
trials  = 200;             % Number of monte carlo trials
snr     = [-100:1:100];     % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
M       = [N/16 N/8 N/4 N/2]; % Filter lengths
a       = 2*fc;               % cyclic frequency two

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*cos(2*pi*fc*nT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pres = zeros(length(snr),length(M));   % probability of resolution

for idx = 1:length(M)   % for each filter length
    for idx2 = 1:length(snr)
    
        success = 0;

        for trial=[1:trials]    % Repeat experiment trials times
    
            disp(['Filter ' int2str(M(idx)) ' SNR ' int2str(snr(idx2)) ' trial ' int2str(trial)]);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add gausian noise to the signals
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            yn = awgn(y,snr(idx2),'measured');
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SCD Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [S1 f1] = scd2(yn,M(idx),'NMVDR',fs,2*fc(1));
            [S2 f2] = scd2(yn,M(idx),'NMVDR',fs,2*fc(2));
            [S3 f3] = scd2(yn,M(idx),'NMVDR',fs,fc(1)+fc(2));

            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq

            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success = success + 1;
            end

        end
    
        pres(idx2,idx)  = success/trials;
    end
end

save 'res_scd03.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation probability of resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,pres(:,1),'o-',snr,pres(:,2),'^-',snr,pres(:,3),'+-',snr,pres(:,4),'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('N/16 length', 'N/8 length', 'N/4 length', 'N/2 length');
title('Probability of resolution of NMVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res03

