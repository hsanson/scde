%
%   File:      powerspectra13.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/6

%   Description:
%        Obtain the variance and vias of MVDR or NMVDR PSD estimators for different SNR
%        and filter lengths.
%
%   Notes:
%        - Tested with Matlab 2007a
%        - Recursive mean and variance calcuation taken from
%               http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Variance.htm
%               http://people.revoledu.com/kardi/tutorial/RecursiveStatistic/Time-Average.htm


clear; clc;
warning off all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = 64;              % Number of observed samples
A       = 1;                % SNR of each sinusoid
fc      = 1.5;              % Frequencies of each sinusoid
fs      = 10;               % Sampling frequency
trials  = 100;                % Number of monte carlo trials
snr     = [-20:5:20];       % SNR in dB
res     = N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
M       = [N/16 N/8 N/4 N/2]; % Filter lengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*cos(2*pi*fc*nT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the theorical PSD of the signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = zeros(N,1);
f = [-N/2:N/2-1]*fs/N;
f1 = [-fliplr(fc); fc];     % frequency bins  
a1 = [fliplr(A) A];         % bins amplitude

for i = 1:length(f1)
    [m,ix] = min(abs(f-f1(i)));          % find index of the closest value to the frequency
    h(ix,:) = (a1(i)^2)/4;               % 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse = zeros(length(snr),length(M));
average = zeros(length(snr),length(M));
variance = zeros(length(snr),length(M));
bias = zeros(length(snr),length(M));

for idx = 1:length(M)   % for each filter length
    for idx2 = 1:length(snr)
    
        mse_sum  = zeros(N,1);
        bias_sum = zeros(N,1);
        mean_sum = zeros(N,1);    
        vari_sum = zeros(N,1);

        for trial=[1:trials]    % Repeat experiment trials times
    
            disp(['Filter ' int2str(M(idx)) ' SNR ' int2str(snr(idx2)) ' trial ' int2str(trial)]);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add gausian noise to the signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            yn = awgn(y,snr(idx2),'measured');
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SCD Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %S = capon(yn,yn,M(idx),'modified');
            S = capon2(yn,yn,M(idx),'modified');
            mse_sum = mse_sum + (S - h).^2;
            bias_sum = bias_sum + (S - h);
            mean_sum = ((trial-1)/trial).*mean_sum + (1/trial).*S;

            if(trial > 1);  
                vari_sum = ((trial-1)/trial).*vari_sum + (1/(trial - 1)).*(S-mean_sum).^2;
            end;
        end
    
        mse(idx2,idx)       = mean(mean(mse_sum./trials)); 
        average(idx2,idx)   = mean(mean(mean_sum)); 
        bias(idx2,idx)      = mean(mean(mean_sum-h)); 
        variance(idx2,idx)  = mean(mean(vari_sum)) + var(var(mean_sum));
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = semilogy(snr,variance(:,1),'o-',snr,variance(:,2),'^-',snr,variance(:,3),'+-',snr,variance(:,4),'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('N/16 length', 'N/8 length', 'N/4 length', 'N/2 length');
title('Variance of NMVDR Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Mean Square Error','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation Bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(2);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = plot(snr,bias(:,1),'o-',snr,bias(:,2),'^-',snr,bias(:,3),'+-',snr,bias(:,4),'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('N/16 length', 'N/8 length', 'N/4 length', 'N/2 length');
title('Bias of NMVDR Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Mean Square Error','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation MSE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(4);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = semilogy(snr,mse(:,1),'o-',snr,mse(:,2),'^-',snr,mse(:,3),'+-',snr,mse(:,4),'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('N/16 length', 'N/8 length', 'N/4 length', 'N/2 length');
title('Mean Squared Error of NMVDR Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Mean Square Error','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis




