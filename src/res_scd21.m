% 
%   File:      res_scd01.m
%   Author(s): Horacio Sanson      
%   Revision : 2007/11/17

%   Description: 
%        Obtain the probability of resolution of WOSA SCD estimators for different frequency separations.
%        
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
N       = 64;               % Number of observed samples
A       = [1 1];            % Amplitude of each sinusoid
fc      = 1;                % Frequencies of each sinusoid
fs      = 10;               % Sampling frequency
trials  = 200;               % Number of monte carlo trials
snr     = [-20:10:20];      % SNR in dB
df      = [0:0.05:1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
M       = N/4;                % Filter lengths
a       = 2*fc;               % cyclic frequency two

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pres = zeros(length(df),length(snr));   % probability of resolution

for idx = 1:length(snr)   % for each filter length
    for idx2 = 1:length(df)
    
        success = 0;

        for trial=[1:trials]    % Repeat experiment trials times
    
            disp(['SNR ' int2str(snr(idx)) ' Freq separation ' num2str(df(idx2)) ' trial ' int2str(trial)]);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate the sinusoidal signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqs = [fc fc+df(idx2)]';
            y = A*cos(2*pi*freqs*nT);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add gausian noise to the signals
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            yn = awgn(y,snr(idx),'measured');
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SCD Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [S1 f1] = scd2(yn,M,'WOSA',fs,2*freqs(1));
            [S2 f2] = scd2(yn,M,'WOSA',fs,2*freqs(2));
            [S3 f3] = scd2(yn,M,'WOSA',fs,freqs(1)+freqs(2));

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

save 'res_scd21.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation probability of resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = plot(df,pres(:,1),'o-',df,pres(:,2),'^-',df,pres(:,3),'+-',df,pres(:,4),'*-',df,pres(:,5),'x-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('SNR -20', 'SNR -10', 'SNR 0', 'SNR 10', 'SNR 20');
title('Probability of resolution of WOSA SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Frequency separation \Deltaf','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res21

