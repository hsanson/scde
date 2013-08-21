%
%   File:      mse02.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/22

%   Description:
%        Ontain the MSE (Minimum Squared Error) of PSD estimators as a function of SNR
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
N   = 64;                     % Number of observed samples
SNR = [-50:5:50];             % Sinusoid range to test
A = 1;                        
f = 0.3;                      % Sinusoid frequency
fs        = 1;                % Sampling frequency
trials    = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts        = 1/fs;                   % Sampling period
nT        = [0:N-1]*ts;             % Time axis

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

mse_wls = zeros(length(SNR),1);
mse_capon1 = zeros(length(SNR),1);
mse_capon2 = zeros(length(SNR),1);
mse_apes = zeros(length(SNR),1);

for idx = 1:length(SNR)   % for each filter length
    
    sum_wls = zeros(N,1);
    sum_capon1 = zeros(N,1);
    sum_capon2 = zeros(N,1);
    sum_apes = zeros(N,1);

    for trial=[1:trials]    % Repeat experiment # trials       

        disp(['SNR ' int2str(SNR(idx)) ' trial ' int2str(trial)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a sinusoidal signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = A*sin(2*pi*f*nT);
        xin = awgn(y,SNR(idx),'measured');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via WLS method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [W,kw] = wls(xin,xin,N/4);
        sum_wls = sum_wls + (W - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via CAPON method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,kc] = capon(xin,xin,N/4);
        sum_capon1 = sum_capon1 + (C - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via CAPON2 method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C2,k2] = capon2(xin,xin,N/8);
        sum_capon2 = sum_capon2 + (C2 - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via APES method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Ca,ka] = apes(xin,xin,N/4);
        sum_apes = sum_apes + (Ca - h).^2;

    end

    mse_wls(idx) = mean(sum_wls./trials); 
    mse_capon1(idx) = mean(sum_capon1./trials); 
    mse_capon2(idx) = mean(sum_capon2./trials); 
    mse_apes(idx) = mean(sum_apes./trials); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation MSE as a function of SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(SNR,mse_wls,'o-',SNR,mse_capon1,'+-',SNR,mse_capon2,'^-',SNR,mse_apes,'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
h = legend('WLS','Capon1','Capon2','APES');
set(h,'FontSize',16);
title('Estimation Mean Squared Error','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('SNR dB','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis

