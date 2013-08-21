%
%   File:      mse_scd06.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/1

%   Description:
%        Ontain the MSE (Minimum Squared Error) of SCD estimators for different SNR
%        and plots them as functions of varying SNR. 
%
%   Notes:
%        - Tested with Matlab 2007a
%        - The APES method is very slow for filter lengths larger than 40
%        - The filter lengths used are based on the best values found empirically
%          in previous simulation (mse_scd01-mse_scd04)
clear; clc;
%warning off all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = 64;               % Number of observed samples
A       = 1;                % SNR of each sinusoid
fc      = 1.5;              % Frequencies of each sinusoid
fs      = 10;               % Sampling frequency
trials  = 200;              % Number of monte carlo trials
snr     = [-20:5:20];       % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;             % Sampling period
nT      = [0:N-1]*ts;       % Time axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*cos(2*pi*fc*nT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the theorical PSD of the noisy signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[h f fa] = scd_theory(A,fc,fs,N,N,0);

h = (A^2)/4;  % the magnitude at the alpha frequency = 2*fc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start monte carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mse = zeros(length(snr),12);
average = zeros(length(snr),12);
variance = zeros(length(snr),12);
bias = zeros(length(snr),12);

for idx2 = 1:length(snr)

    mse_sum_wosa  =0;
    bias_sum_wosa =0;
    mean_sum_wosa =0;    
    vari_sum_wosa =0;
    mse_sum_wosa2  =0;
    bias_sum_wosa2 =0;
    mean_sum_wosa2 =0;    
    vari_sum_wosa2 =0;
    mse_sum_wosa3  =0;
    bias_sum_wosa3 =0;
    mean_sum_wosa3 =0;    
    vari_sum_wosa3 =0;

    mse_sum_mvdr  =0;
    bias_sum_mvdr =0;
    mean_sum_mvdr =0;    
    vari_sum_mvdr =0;
    mse_sum_mvdr2  =0;
    bias_sum_mvdr2 =0;
    mean_sum_mvdr2 =0;    
    vari_sum_mvdr2 =0;
    mse_sum_mvdr3  =0;
    bias_sum_mvdr3 =0;
    mean_sum_mvdr3 =0;    
    vari_sum_mvdr3 =0;

    mse_sum_nmvdr  =0;
    bias_sum_nmvdr =0;
    mean_sum_nmvdr =0;    
    vari_sum_nmvdr =0;
    mse_sum_nmvdr2  =0;
    bias_sum_nmvdr2 =0;
    mean_sum_nmvdr2 =0;    
    vari_sum_nmvdr2 =0;
    mse_sum_nmvdr3  =0;
    bias_sum_nmvdr3 =0;
    mean_sum_nmvdr3 =0;    
    vari_sum_nmvdr3 =0;
    
    mse_sum_apes  =0;
    bias_sum_apes =0;
    mean_sum_apes =0;    
    vari_sum_apes =0;
    mse_sum_apes2  =0;
    bias_sum_apes2 =0;
    mean_sum_apes2 =0;    
    vari_sum_apes2 =0;
    mse_sum_apes3  =0;
    bias_sum_apes3 =0;
    mean_sum_apes3 =0;    
    vari_sum_apes3 =0;

    for trial=[1:trials]    % Repeat experiment trials times

        disp(['SNR ' int2str(snr(idx2)) ' trial ' int2str(trial)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add gausian noise to the signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xin = awgn(y,snr(idx2),'measured');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SCD Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C f] = scd2(xin,N,'WOSA',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_wosa = mse_sum_wosa + (S - h).^2;
        bias_sum_wosa = bias_sum_wosa + (S - h);
        mean_sum_wosa = ((trial-1)/trial).*mean_sum_wosa + (1/trial).*S;
        if(trial > 1);  
            vari_sum_wosa = ((trial-1)/trial).*vari_sum_wosa + (1/(trial - 1)).*(S-mean_sum_wosa).^2;
        end;

        [C f] = scd2(xin,N/2,'WOSA',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_wosa2 = mse_sum_wosa2 + (S - h).^2;
        bias_sum_wosa2 = bias_sum_wosa2 + (S - h);
        mean_sum_wosa2 = ((trial-1)/trial).*mean_sum_wosa2 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_wosa2 = ((trial-1)/trial).*vari_sum_wosa2 + (1/(trial - 1)).*(S-mean_sum_wosa2).^2;
        end;

        [C f] = scd2(xin,N/4,'WOSA',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_wosa3 = mse_sum_wosa3 + (S - h).^2;
        bias_sum_wosa3 = bias_sum_wosa3 + (S - h);
        mean_sum_wosa3 = ((trial-1)/trial).*mean_sum_wosa3 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_wosa3 = ((trial-1)/trial).*vari_sum_wosa3 + (1/(trial - 1)).*(S-mean_sum_wosa3).^2;
        end;

        [C f] = scd2(xin,N/2,'MVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_mvdr = mse_sum_mvdr + (S - h).^2;
        bias_sum_mvdr = bias_sum_mvdr + (S - h);
        mean_sum_mvdr = ((trial-1)/trial).*mean_sum_mvdr + (1/trial).*S;
        if(trial > 1);  
            vari_sum_mvdr = ((trial-1)/trial).*vari_sum_mvdr + (1/(trial - 1)).*(S-mean_sum_mvdr).^2;
        end;

        [C f] = scd2(xin,N/4,'MVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_mvdr2 = mse_sum_mvdr2 + (S - h).^2;
        bias_sum_mvdr2 = bias_sum_mvdr2 + (S - h);
        mean_sum_mvdr2 = ((trial-1)/trial).*mean_sum_mvdr2 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_mvdr2 = ((trial-1)/trial).*vari_sum_mvdr2 + (1/(trial - 1)).*(S-mean_sum_mvdr2).^2;
        end;

        [C f] = scd2(xin,N/8,'MVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                    % find index of the closest zero freq
        S = real(C(k));

        mse_sum_mvdr3 = mse_sum_mvdr3 + (S - h).^2;
        bias_sum_mvdr3 = bias_sum_mvdr3 + (S - h);
        mean_sum_mvdr3 = ((trial-1)/trial).*mean_sum_mvdr3 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_mvdr3 = ((trial-1)/trial).*vari_sum_mvdr3 + (1/(trial - 1)).*(S-mean_sum_mvdr3).^2;
        end;


        [C f] = scd2(xin,N/2,'NMVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_nmvdr = mse_sum_nmvdr + (S - h).^2;
        bias_sum_nmvdr = bias_sum_nmvdr + (S - h);
        mean_sum_nmvdr = ((trial-1)/trial).*mean_sum_nmvdr + (1/trial).*S;
        if(trial > 1);  
            vari_sum_nmvdr = ((trial-1)/trial).*vari_sum_nmvdr + (1/(trial - 1)).*(S-mean_sum_nmvdr).^2;
        end;

        [C f] = scd2(xin,N/4,'NMVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_nmvdr2 = mse_sum_nmvdr2 + (S - h).^2;
        bias_sum_nmvdr2 = bias_sum_nmvdr2 + (S - h);
        mean_sum_nmvdr2 = ((trial-1)/trial).*mean_sum_nmvdr2 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_nmvdr2 = ((trial-1)/trial).*vari_sum_nmvdr2 + (1/(trial - 1)).*(S-mean_sum_nmvdr2).^2;
        end;

        [C f] = scd2(xin,N/16,'NMVDR',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_nmvdr3 = mse_sum_nmvdr3 + (S - h).^2;
        bias_sum_nmvdr3 = bias_sum_nmvdr3 + (S - h);
        mean_sum_nmvdr3 = ((trial-1)/trial).*mean_sum_nmvdr3 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_nmvdr3 = ((trial-1)/trial).*vari_sum_nmvdr3 + (1/(trial - 1)).*(S-mean_sum_nmvdr3).^2;
        end;


        [C f] = scd2(xin,N/2,'APES',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_apes = mse_sum_apes + (S - h).^2;
        bias_sum_apes = bias_sum_apes + (S - h);
        mean_sum_apes = ((trial-1)/trial).*mean_sum_apes + (1/trial).*S;
        if(trial > 1);  
            vari_sum_apes = ((trial-1)/trial).*vari_sum_apes + (1/(trial - 1)).*(S-mean_sum_apes).^2;
        end;

        [C f] = scd2(xin,N/4,'APES',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_apes2 = mse_sum_apes2 + (S - h).^2;
        bias_sum_apes2 = bias_sum_apes2 + (S - h);
        mean_sum_apes2 = ((trial-1)/trial).*mean_sum_apes2 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_apes2 = ((trial-1)/trial).*vari_sum_apes2 + (1/(trial - 1)).*(S-mean_sum_apes2).^2;
        end;

        [C f] = scd2(xin,N/8,'APES',fs,2*fc);   % Obtain the SCD at alpha = 2*fc
        [m k] = min(abs(f));                 % find index of the closest zero freq
        S = real(C(k));

        mse_sum_apes3 = mse_sum_apes3 + (S - h).^2;
        bias_sum_apes3 = bias_sum_apes3 + (S - h);
        mean_sum_apes3 = ((trial-1)/trial).*mean_sum_apes3 + (1/trial).*S;
        if(trial > 1);  
            vari_sum_apes3 = ((trial-1)/trial).*vari_sum_apes3 + (1/(trial - 1)).*(S-mean_sum_apes3).^2;
        end;
    end

    mse(idx2,1)       = mse_sum_wosa./trials; 
    average(idx2,1)   = mean_sum_wosa; 
    bias(idx2,1)      = bias_sum_wosa./trials; 
    variance(idx2,1)  = vari_sum_wosa ;

    mse(idx2,2)       = mse_sum_wosa2./trials; 
    average(idx2,2)   = mean_sum_wosa2; 
    bias(idx2,2)      = bias_sum_wosa2./trials; 
    variance(idx2,2)  = vari_sum_wosa2 ;

    mse(idx2,3)       = mse_sum_wosa3./trials; 
    average(idx2,3)   = mean_sum_wosa3; 
    bias(idx2,3)      = bias_sum_wosa3./trials; 
    variance(idx2,3)  = vari_sum_wosa3;

    mse(idx2,4)       = mse_sum_mvdr./trials; 
    average(idx2,4)   = mean_sum_mvdr; 
    bias(idx2,4)      = bias_sum_mvdr./trials; 
    variance(idx2,4)  = vari_sum_mvdr ;

    mse(idx2,5)       = mse_sum_mvdr2./trials; 
    average(idx2,5)   = mean_sum_mvdr2; 
    bias(idx2,5)      = bias_sum_mvdr2./trials; 
    variance(idx2,5)  = vari_sum_mvdr2; 
    
    mse(idx2,6)       = mse_sum_mvdr3./trials; 
    average(idx2,6)   = mean_sum_mvdr3; 
    bias(idx2,6)      = bias_sum_mvdr3./trials; 
    variance(idx2,6)  = vari_sum_mvdr3;
    
    mse(idx2,7)       = mse_sum_nmvdr./trials ;
    average(idx2,7)   = mean_sum_nmvdr ;
    bias(idx2,7)      = bias_sum_nmvdr./trials; 
    variance(idx2,7)  = vari_sum_nmvdr;
 
    mse(idx2,8)       = mse_sum_nmvdr2./trials ;
    average(idx2,8)   = mean_sum_nmvdr2 ;
    bias(idx2,8)      = bias_sum_nmvdr2./trials;
    variance(idx2,8)  = vari_sum_nmvdr2;

    mse(idx2,9)       = mse_sum_nmvdr3./trials ;
    average(idx2,9)   = mean_sum_nmvdr3 ;
    bias(idx2,9)      = bias_sum_nmvdr3./trials;
    variance(idx2,9)  = vari_sum_nmvdr3;

    mse(idx2,10)       = mse_sum_apes./trials; 
    average(idx2,10)   = mean_sum_apes; 
    bias(idx2,10)      = bias_sum_apes./trials; 
    variance(idx2,10)  = vari_sum_apes;

    mse(idx2,11)       = mse_sum_apes2./trials; 
    average(idx2,11)   = mean_sum_apes2; 
    bias(idx2,11)      = bias_sum_apes2./trials; 
    variance(idx2,11)  = vari_sum_apes2;

    mse(idx2,12)       = mse_sum_apes3./trials; 
    average(idx2,12)   = mean_sum_apes3; 
    bias(idx2,12)      = bias_sum_apes3./trials; 
    variance(idx2,12)  = vari_sum_apes3;
end

save 'mse_scd06.mat'

%bias = bias.^2;    % Square the bias to remove negative values... if not we get gaps in the logarithmic graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(1);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,1),'o-',snr,variance(:,2),'^-',snr,variance(:,3),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N', 'WOSA N/2', 'WOSA N/4');
title('Variance of the WOSA SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance01

fh = figure(2);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,4),'o-',snr,variance(:,5),'^-',snr,variance(:,6),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('MVDR N/2', 'MVDR N/4', 'MVDR N/8');
title('Variance of the MVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance02

fh = figure(3);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,7),'o-',snr,variance(:,8),'^-',snr,variance(:,9),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8');
title('Variance of the NMVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance03

fh = figure(4);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,10),'o-',snr,variance(:,11),'^-',snr,variance(:,12),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('APES N/2', 'APES N/4', 'APES N/8');
title('Variance of the APES SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation Bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(5);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,1),'o-',snr,bias(:,2),'^-',snr,bias(:,3),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N', 'WOSA N/2', 'WOSA N/4');
title('Bias of the WOSA SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias01

fh = figure(6);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,4),'o-',snr,bias(:,5),'^-',snr,bias(:,6),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('MVDR N/2', 'MVDR N/4', 'MVDR N/8');
title('Bias of the MVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias02

fh = figure(7);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,7),'o-',snr,bias(:,8),'^-',snr,bias(:,9),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8');
title('Bias of the NMVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias03

fh = figure(8);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,10),'o-',snr,bias(:,11),'^-',snr,bias(:,12),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('APES N/2', 'APES N/4', 'APES N/8');
title('Bias of the APES SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot them all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure(9);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,1),'o-',snr,variance(:,2),'o--',snr,variance(:,3),'o-.',snr,variance(:,4),'s-',snr,variance(:,5),'s--',snr,variance(:,6),'s-.',snr,variance(:,7),'^-',snr,variance(:,8),'^--',snr,variance(:,9),'^-.',snr,variance(:,10),'h-',snr,variance(:,11),'h--',snr,variance(:,12),'h-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8','NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8','APES N/2', 'APES N/4', 'APES N/8');
title('Variance Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance06

fh = figure(10);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,1),'o-',snr,bias(:,2),'o--',snr,bias(:,3),'o-.',snr,bias(:,4),'s-',snr,bias(:,5),'s--',snr,bias(:,6),'s-.',snr,bias(:,7),'^-',snr,bias(:,8),'^--',snr,bias(:,9),'^-.',snr,bias(:,10),'h-',snr,bias(:,11),'h--',snr,bias(:,12),'h-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8','NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8','APES N/2', 'APES N/4', 'APES N/8');
title('Bias Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias06

fh = figure(11);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,variance(:,1),'o-',snr,variance(:,2),'o--',snr,variance(:,3),'o-.',snr,variance(:,4),'s-',snr,variance(:,5),'s--',snr,variance(:,6),'s-.',snr,variance(:,7),'^-',snr,variance(:,8),'^--',snr,variance(:,9),'^-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8','NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8');
title('Variance Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('Variance','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_variance05

fh = figure(12);
set(fh, 'color', 'white'); % sets the color to white
fig1 = plot(snr,bias(:,1),'o-',snr,bias(:,2),'o--',snr,bias(:,3),'o-.',snr,bias(:,4),'s-',snr,bias(:,5),'s--',snr,bias(:,6),'s-.',snr,bias(:,8),'^--',snr,bias(:,9),'^-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8', 'NMVDR N/4', 'NMVDR N/8');
title('Bias Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('Bias','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_bias05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimation MSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(13);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,1),'o-',snr,mse(:,2),'^-',snr,mse(:,3),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('WOSA N', 'WOSA N/2', 'WOSA N/4');
title('MSE of the WOSA SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse01

fh = figure(14);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,4),'o-',snr,mse(:,5),'^-',snr,mse(:,6),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('MVDR N/2', 'MVDR N/4', 'MVDR N/8');
title('MSE of the MVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse02

fh = figure(15);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,7),'o-',snr,mse(:,8),'^-',snr,mse(:,9),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8');
title('MSE of the NMVDR SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse03

fh = figure(16);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,10),'o-',snr,mse(:,11),'^-',snr,mse(:,12),'s-');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('APES N/2', 'APES N/4', 'APES N/8');
title('MSE of the APES SCD Estimator','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse04


fh = figure(18);
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,1),'o-',snr,mse(:,2),'o--',snr,mse(:,3),'o-.',snr,mse(:,4),'s-',snr,mse(:,5),'s--',snr,mse(:,6),'s-.',snr,mse(:,7),'^-',snr,mse(:,8),'^--',snr,mse(:,9),'^-.',snr,variance(:,10),'h-',snr,variance(:,11),'h--',snr,variance(:,12),'h-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8','NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8','APES N/2', 'APES N/4', 'APES N/8');
title('MSE Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse06

fh = figure(19);
%set(0,'DefaultAxesColorOrder',[0 0 0]);    % makes the graph black and white
set(fh, 'color', 'white'); % sets the color to white
fig1 = semilogy(snr,mse(:,1),'o-',snr,mse(:,2),'o--',snr,mse(:,3),'o-.',snr,mse(:,4),'s-',snr,mse(:,5),'s--',snr,mse(:,6),'s-.',snr,mse(:,7),'^-',snr,mse(:,8),'^--',snr,mse(:,9),'^-.');
set(fig1, 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend( 'WOSA N', 'WOSA N/2', 'WOSA N/4','MVDR N/2', 'MVDR N/4', 'MVDR N/8','NMVDR N/2', 'NMVDR N/4', 'NMVDR N/8');
title('MSE Comparison of SCD Estimators','FontSize',16,'FontWeight', 'bold');
ylabel('MSE','FontSize',16,'FontWeight', 'bold');
xlabel('Signal to Noise Ratio (dB)','FontSize',16,'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis
print -deps -painters  img/scd_mse05
