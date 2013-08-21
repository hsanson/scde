% 
%   File:      res_scd28.m
%   Author(s): Horacio Sanson      
%   Revision : 2007/11/17

%   Description: 
%        Obtain the probability of resolution of SCD estimators for different frequency separations.
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
trials  = 200;              % Number of monte carlo trials
snr     = 0;                % SNR in dB
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

pres = zeros(length(df),8);   % probability of resolution

for idx = 1:length(df)   % for each filter length

        success_wosa = 0;
        success_wosa2 = 0;
        success_mvdr = 0;
        success_mvdr2 = 0;
        success_nmvdr = 0;
        success_nmvdr2 = 0;
        success_apes = 0;
        success_apes2 = 0;

        for trial=[1:trials]    % Repeat experiment trials times
    
            disp(['Freq separation ' num2str(df(idx)) ' trial ' int2str(trial)]);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate the sinusoidal signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqs = [fc fc+df(idx)]';
            y = A*cos(2*pi*freqs*nT);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add gausian noise to the signals
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            yn = awgn(y,snr,'measured');
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SCD Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [S1 f1] = scd2(yn,N/4,'WOSA',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/4,'WOSA',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/4,'WOSA',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_wosa = success_wosa + 1;
            end
    
            [S1 f1] = scd2(yn,N,'WOSA',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N,'WOSA',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N,'WOSA',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_wosa2 = success_wosa2 + 1;
            end

            [S1 f1] = scd2(yn,N/4,'MVDR',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/4,'MVDR',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/4,'MVDR',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_mvdr = success_mvdr + 1;
            end
    
            [S1 f1] = scd2(yn,N/8,'MVDR',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/8,'MVDR',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/8,'MVDR',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_mvdr2 = success_mvdr2 + 1;
            end

            [S1 f1] = scd2(yn,N/4,'NMVDR',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/4,'NMVDR',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/4,'NMVDR',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_nmvdr = success_nmvdr + 1;
            end
    
            [S1 f1] = scd2(yn,N/8,'NMVDR',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/8,'NMVDR',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/8,'NMVDR',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_nmvdr2 = success_nmvdr2 + 1;
            end
    
            [S1 f1] = scd2(yn,N/4,'APES',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/4,'APES',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/4,'APES',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_apes = success_apes + 1;
            end

            [S1 f1] = scd2(yn,N/8,'APES',fs,2*freqs(1));
            [S2 f2] = scd2(yn,N/8,'APES',fs,2*freqs(2));
            [S3 f3] = scd2(yn,N/8,'APES',fs,freqs(1)+freqs(2));
    
            % Find the bin at frequency 0
            [mag1 k1] = min(abs(f1));                 % find index of the closest zero freq
            [mag2 k2] = min(abs(f2));                 % find index of the closest zero freq
            [mag3 k3] = min(abs(f3));                 % find index of the closest zero freq
    
            if ((1/2)*(abs(S2(k2))+abs(S1(k1))) - abs(S3(k3))) > 0
                success_apes2 = success_apes2 + 1;
            end
        end
    
        pres(idx,1)  = success_wosa/trials;
        pres(idx,2)  = success_wosa2/trials;
        pres(idx,3)  = success_mvdr/trials;
        pres(idx,4)  = success_mvdr2/trials;
        pres(idx,5)  = success_nmvdr/trials;
        pres(idx,6)  = success_nmvdr2/trials;
        pres(idx,7)  = success_apes/trials;  
        pres(idx,8)  = success_apes2/trials;
     
end

save 'res_scd30.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation probability of resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = plot(df,pres(:,1),'o-',df,pres(:,2),'^-',df,pres(:,3),'+-',df,pres(:,4),'*-',df,pres(:,5),'x-',df,pres(:,6),'s-',df,pres(:,7),'d-',df,pres(:,8),'v-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N/4', 'WOSA N', 'MVDR N/4', 'MVDR N/8', 'NMVDR N/4', 'NMVDR N/8', 'APES N/4', 'APES N/8');
title('Probability of resolution of SCD Estimators (0dB SNR)','FontSize',16,'FontWeight', 'bold');
ylabel('Probability of Resolution P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Cyclic Frequency Separation \partial\alpha','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res29


fh = figure(2);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = plot(df,pres(:,1),'o-',df,pres(:,2),'^-',df,pres(:,3),'+-',df,pres(:,4),'*-',df,pres(:,5),'x-',df,pres(:,6),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N/4', 'WOSA N', 'MVDR N/4', 'MVDR N/8', 'NMVDR N/4', 'NMVDR N/8');
title('Probability of resolution of SCD Estimators (0dB SNR)','FontSize',16,'FontWeight', 'bold');
ylabel('Probability of Resolution P(\Gamma > 0)','FontSize',16,'FontWeight', 'bold');
xlabel('Cyclic Frequency Separation \partial\alpha','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_res30


