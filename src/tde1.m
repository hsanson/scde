%
%   File:      tde1.m
%   Author(s): Horacio Sanson
%
%   Description:
%     Time Delay estimation for OFDM signal.
%     Uses a basic 1014-FFT OFDMA-PHY tranceiver

clear; clc;

C = 299792458;		% Speed of Ligth constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDMA-PHY Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Primitive Parameters (8.4.2.3)
Bw = 5e6;                                   % Nominal Bandwidth 
n = 8/7;                                    % Oversampling Factor
G = 1/8;                                    % CP time to Tb factor
Nfft = 1024;                                % FFT Size

% Derived Parameters (8.4.2.4)
Fs = floor(n*Bw/8000)*8000;                 % Sampling Freq
deltaF = Fs/Nfft;                           % Subcarrier Spacing
Tb = 1/deltaF;                              % Useful symbol time
Tg = G*Tb;                                  % Guard Interval
Ts = Tb + Tg;                               % OFDMA Symbol time
T = Tb/Nfft;                                % Sampling Time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Modulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 4;                          % Alphabet Size (QPSK)
k = log2(M);                    % Number of bits per symbol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pilot/Data Subcarrier Allocation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation for Uplink PUSC 1024-FFT (Table 313a)
Nleft = 92;                     % Left guard subcarriers 
Nright = 91;                    % Rigth guard subcarriers
Nused = 841;                    % Used subcarriers (Ndata+Ndc+Npilot)
Ndc = 1;                        % DC subcarrier, index 512
Ndata = 840;                    % Assign all data subcarriers (to simplify)
Npilot = 0;                     % Num pilot subcarriers

UsedSubcIndex = [93:511 513:933]';
PilotSubcIndex = []';			% We won't use pilots in this simulation	

%%%%%%%%%%%%%%%%%%%%
% SIMULATION START %
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ofdm_symbols = 1;                        % Number of OFDM symbols to transmit
packet_length = n_ofdm_symbols*k*Ndata;     % Num de bits per packet

x = randint(1,packet_length);               % Generate random data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QPSK Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbit = reshape(x,k,length(x)/k)';
xsym = bi2de(xbit, 'left-msb');
xmod = pskmod(xsym',M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate Subcarriers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mod_syms = size(xmod,2);
n_ofdm_syms = n_mod_syms/Ndata;

% Create a new matrix and map xmod based on the UsedSubcIndex
xofdm = zeros(Nfft, n_ofdm_syms);
xofdm(UsedSubcIndex,:) = reshape(xmod, Ndata, n_ofdm_syms);

% Fill the Pilot Subcarriers with the pattern (not implemented)
xofdm(PilotSubcIndex,:) = -1;   % Just fill them with -1's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to Time Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtime = ifft(xofdm,Nfft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the Cyclic Prefix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_size = Nfft*G;
tail_start = Nfft - cp_size + 1;
tail_end = Nfft;

xcptime = [ xtime(tail_start:tail_end,:) ; xtime ];  % Add cyclic CP
%xcptime = [zeros(cp_size,n_ofdm_symbols); xtime];   % Add zeros CP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P/S conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrf = xcptime(:).';     % Note: remember that .' and ' are different. Using
                        % only ' will conjugate complex values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create delayed version of the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = 60;		      % samples delay

yrf = filter([1 1], 1, xrf);  % Create delayed version of xrf
yrf = [zeros(1,delay) yrf];
yrf = yrf(1:length(xrf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit both signals throught an AWGN Channel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ebno = -5;
snr = ebno + 10*log10(k);


%delays = [Ts/150 Ts/20];
%powers = [0 -3.2];
%rayleigh = rayleighchan(1/Fs,40,delays,powers);
%xray = filter(rayleigh,xrf);
%yray = filter(rayleigh,yrf);


% Add AWGN noise to the xrf signal
xnoisy = awgn(xrf, snr, 'measured');
% Add AWGN noise to the yrf signal
ynoisy = awgn(yrf, snr, 'measured');

%xnoisy=xrf;
%ynoisy=yrf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the spectra.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Sxx, Fxx] = pwelch(xnoisy,[],[],Nfft, Fs);
[Syy, Fyy] = pwelch(ynoisy,[],[],Nfft, Fs);
[Sxy, Fxy] = cpsd(xnoisy,ynoisy,[],[], Nfft, Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the Generalized Correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[corr1, t1] = GCC('unfiltered',Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);
[corr2, t2] = GCC('scot', Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);
[corr3, t3] = GCC('roth', Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);
[corr4, t4] = GCC('phat', Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);
[corr5, t5] = GCC('cps-m', Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);
[corr6, t6] = GCC('ht', Sxx, Syy, Sxy, Fs, length(xnoisy), Nfft);

figure(1);
subplot(321);
plot(t1,corr1);
%axis([-1 1 -1 1]);
[m,k] = max(corr1);
k = t1(k);
title('Unfiltered GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

subplot(322);
plot(t2,corr2);
%axis([-30e-6 30e-6 -1.5 1.5]);
[m,k] = max(corr2);
k = t2(k);
title('SCOT GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

subplot(323);
plot(t3,corr3);
%axis([-1 1 -1 1]);
[m,k] = max(corr3);
k = t3(k);
title('ROTH GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

subplot(324);
plot(t4,corr4);
%axis([-1 1 -1 1]);
[m,k] = max(corr4);
k = t4(k);
title('PHAT GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

subplot(325);
plot(t5,corr5);
%axis([-1 1 -1 1]);
[m,k] = max(corr5);
k = t5(k);
title('cps-m GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

subplot(326);
plot(t6,corr6);
%axis([-1 1 -1 1]);
[m,k] = max(corr6);
k = t6(k);
title('HT GCC');
text(1e-4,0.5,strcat(num2str(k),'sec'));

% Plot power spectrum of OFDM signal
%figure(2);
%subplot(211);
%[P,F] = pwelch(xrf,[],[],Nfft,Fs);
%semilogy(F,P);
%title('Tx OFDM Spectrum');
%subplot(212);
%[P,F] = pwelch(xnoisy,[],[],Nfft,Fs);
%semilogy(F,P);
%title('Tx OFDM Spectrum with AWGN');


