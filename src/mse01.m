%
%   File:      mse01.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/22

%   Description:
%        Obtain the MSE (Minimum Squared Error) of PSD estimators for different filter lengths
%
%   Notes:
%        - Tested with Matlab 7.0.0.19901 (R14)
%        - The APES method is very slow for filter lengths larger than 40

clear; clc;
warning off all;  % Needed to supress some singular matrix warnings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N   = 64;                      % Number of observed samples
M   = [8:4:N/2];               % Filter lengths to test
A   = [1 1];                   % SNR of each sinusoid
f = [0.3;0.2];                 % Frequencies of each sinusoid
fs        = 1;                 % Sampling frequency
trials    = 10;
snr       = 35;                % SNR in dB

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

mse_wls = zeros(length(M),1);
mse_capon1 = zeros(length(M),1);
mse_capon2 = zeros(length(M),1);
mse_apes = zeros(length(M),1);

for idx = 1:length(M)   % for each filter length
    
    sum_wls = zeros(N,1);
    sum_capon1 = zeros(N,1);
    sum_capon2 = zeros(N,1);
    sum_apes = zeros(N,1);

    for trial=[1:trials]    % Repeat experiment # trials       

        disp(['Filter ' int2str(M(idx)) ' trial ' int2str(trial)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add gausian noise to the signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xin = awgn(y,snr,'measured');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via WLS method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [W,kw] = wls(xin,xin,M(idx));
        sum_wls = sum_wls + (W - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via CAPON method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,kc] = capon(xin,xin,M(idx));
        sum_capon1 = sum_capon1 + (C - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via CAPON2 method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C2,k2] = capon2(xin,xin,M(idx));
        sum_capon2 = sum_capon2 + (C2 - h).^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation via APES method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Ca,ka] = apes(xin,xin,M(idx));
        sum_apes = sum_apes + (Ca - h).^2;

    end

    mse_wls(idx) = mean(sum_wls./trials); 
    mse_capon1(idx) = mean(sum_capon1./trials); 
    mse_capon2(idx) = mean(sum_capon2./trials); 
    mse_apes(idx) = mean(sum_apes./trials); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation MSE as a function of filter length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(M,mse_wls,'o-',M,mse_capon1,'+-',M,mse_capon2,'^-',M,mse_apes,'*-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
h = legend('WOSA (FFT)','Capon1 (MLM)','Capon2 (NMLM)','APES');
set(h,'FontSize',16);
title('Estimated PSD Mean Squared Error','FontSize',16,'FontWeight', 'bold');
ylabel('Mean Square Error','FontSize',16,'FontWeight', 'bold');
xlabel('Filter Length','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the PSD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
grid on;
plot(fh,h,'--',kw*fs,W,kc*fs,C,k2*fs,C2,ka*fs,Ca);
legend('Theorical','WLS','Capon','Capon2','APES');
title('Spectral Density');
ylabel('PSD');
xlabel('Hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the PSD in dB scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
grid on;
plot(fh,hdb,'--',kw*fs,10*log10(W),kc*fs,10*log10(C),k2*fs,10*log10(C2),ka*fs,10*log10(Ca));
legend('Theorical','WLS','Capon','Capon2','APES');
title('Spectral Density');
ylabel('PSD db');
xlabel('Hz');

