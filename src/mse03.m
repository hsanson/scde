%
%   File:      mse03.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/22

%   Description:
%        Obtain the MSE (Minimum Squared Error) of APES PSD estimators for different filter lengths
%        and plots them as functions of varying SNR
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14)
%        - The APES method is very slow for filter lengths larger than 40

clear; clc;
warning off all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N   = 64;                      % Number of observed samples
M   = [4:4:N/2];               % Filter lengths to test
A   = [1 1];                   % SNR of each sinusoid
f = [0.3;0.1];                 % Frequencies of each sinusoid
fs        = 1;                 % Sampling frequency
trials    = 10;
snr       = [-20:5:50];        % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts        = 1/fs;                   % Sampling period
nT        = [0:N-1]*ts;             % Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*sin(2*pi*f*nT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the theorical PSD of the generated signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h   = zeros(1,N);
fh  = [-N/2:N/2-1]*fs/N;
f1 = [-fliplr(f); f];
a1 = [fliplr(A) A];
for i = 1:length(f1)
    [m,ix] = min(abs(fh-f1(i)));   % find index of the closest value to the frequency
    h(ix) = (a1(i)^2)/2;           % theorical PSD
end
h = h';
hdb = 10*log10(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mse = zeros(length(snr),length(M));

for idx = 1:length(M)   % for each filter length
   
    for idx2 = 1:length(snr)

        mse_sum = zeros(N,1);

        for trial=[1:trials]    % Repeat experiment # trials       

            disp(['Filter ' int2str(M(idx)) ' SNR ' int2str(snr(idx2)) ' trial ' int2str(trial)]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add gausian noise to the signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xin = awgn(y,snr(idx2),'measured');


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimation via CAPON method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [C,kc] = apes(xin,xin,M(idx));
            mse_sum = mse_sum + (C - h).^2;

        end

        mse(idx2,idx) = mean(mse_sum./trials); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation MSE as a function of filter length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
fig1 = semilogy(snr,mse);
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend(int2str(M(:)));
title('Estimated PSD Mean Squared Error','FontSize',16,'FontWeight', 'bold');
ylabel('Mean Square Error','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis


