%
%   File:      powerspectra10.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17

%   Description:
%        Demo showing how to estimate the Spectral Correlation Density of a sinusoid signal
%
%   Notes:
%       - Tested with Matlab 7.0.0.19901 (R14) and Octave 2.9.9
%       - Tested with Matlab 2007a

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 1;                         % Signal freq
fs = 10;                        % Sampling Freq
N = 64;                        % Number of signal samples
snr = 10;                       % Signal to noise ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 1/fs;                      % Sampling period
n = [0:N-1]'*Ts;                % Time vector
fftlength = N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn = cos(2*pi*fc.*n);
sn_noisy = awgn(sn,snr,'measured');  % Add AWGN noise

% Theorical SCD of AM signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the theorical SCD of the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S = zeros(length(f),length(a));
%for i=1:length(a)
%    Q1 = sin(pi*td*(f + fc + a(i)/2))./(pi*(f + fc + a(i)/2));
%    Q2 = sin(pi*td*(f + fc - a(i)/2))./(pi*(f + fc - a(i)/2));
%    Q3 = sin(pi*td*(f - fc + a(i)/2))./(pi*(f - fc + a(i)/2));
%    Q4 = sin(pi*td*(f - fc - a(i)/2))./(pi*(f - fc - a(i)/2));
%
%    S1(:,i) =  (1/(2*td)) * ( (Q1.*Q2 + Q3.*Q4) + (Q1.*Q4) + (Q2.*Q3) );
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the SCD using different methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h f fa] = scd_theory(1,fc,fs,N,N,0);
[S2 f a] = scd(sn_noisy,N/2,'WOSA',fs);
[C1 k1 a1] = scd(sn_noisy,N/2,'MVDR',fs); 
[C2 k2 a2] = scd(sn_noisy,N/4,'NMVDR',fs); 
[Ca ka aa] = scd(sn_noisy,N/2,'APES',fs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot only the Capon1 method only for publication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure(2);
set(fh, 'color', 'white'); % sets the color to white

subplot(221);
mesh(a,f,abs(h));
title('Theorical SCD','FontSize',16,'FontWeight', 'bold');
zlabel('|PSD|','FontSize',12,'FontWeight', 'bold');
xlabel('\alpha (Hz)','FontSize',12,'FontWeight', 'bold');
ylabel('freq (Hz)','FontSize',12,'FontWeight', 'bold');
subplot(222);
mesh(a,f,abs(S2));
title('WOSA SCD','FontSize',16,'FontWeight', 'bold');
zlabel('|PSD|','FontSize',12,'FontWeight', 'bold');
xlabel('\alpha (Hz)','FontSize',12,'FontWeight', 'bold');
ylabel('freq (Hz)','FontSize',12,'FontWeight', 'bold');
subplot(223);
mesh(a,f,abs(C1));
title('MVDR SCD','FontSize',16,'FontWeight', 'bold');
zlabel('|PSD|','FontSize',12,'FontWeight', 'bold');
xlabel('\alpha (Hz)','FontSize',12,'FontWeight', 'bold');
ylabel('freq (Hz)','FontSize',12,'FontWeight', 'bold');
subplot(224);
mesh(a,f,abs(C2));
title('NMVDR SCD','FontSize',16,'FontWeight', 'bold');
zlabel('|PSD|','FontSize',12,'FontWeight', 'bold');
xlabel('\alpha (Hz)','FontSize',12,'FontWeight', 'bold');
ylabel('freq (Hz)','FontSize',12,'FontWeight', 'bold');

print -deps -painters  img/cyclicspec

figure(3);

subplot(2,2,1);
contour(aa,ka,abs(Ca));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('Cyclic Spectrum (APES)','FontSize', 16);

subplot(2,2,2);
contour(a,f,abs(S2));
title('Cyclic Spectrum (WLS)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(2,2,3);
contour(a1,k1,abs(C1));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('Cyclic Spectrum (CAPON-1)','FontSize', 16);

subplot(2,2,4);
contour(a2,k2,abs(C2));
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
title('Cyclic Spectrum (CAPON-2)','FontSize', 16);

figure(4);

subplot(2,2,1);
mesh(aa,ka,abs(Ca));
title('Cyclic Spectrum (APES)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(2,2,2);
mesh(a,f,abs(S2));
title('Cyclic Spectrum (Welch)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(2,2,3);
mesh(a1,k1,abs(C1));
title('Cyclic Spectrum (CAPON-1)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

subplot(2,2,4);
mesh(a2,k2,abs(C2));
title('Cyclic Spectrum (CAPON-2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

%set(gca, 'Box', 'off','TickDir', 'out', 'FontSize',16 ); % here gca means get current axis

