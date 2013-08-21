%
%   File:      ofdm_test4.m
%   Author(s): Horacio Sanson

%   Description:
%     Simulate a basic 1014-FFT OFDMA-PHY tranceiver and calculate it's BER
%     under several Eb/No values.
%
%     The BER simulation of a baseband OFDM system using QPSK should coincide
%     with the BER simulation of single carrier QPSK.
%
%     With this simulation we validate our model.

clear; clc;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BER Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbNoVec = [0:0.5:10];
number_of_errors = zeros(1, length(EbNoVec));
bit_error_rate = zeros(1, length(EbNoVec));

%%%%%%%%%%%%%%%%%%%%
% SIMULATION START %
%%%%%%%%%%%%%%%%%%%%

for idxEbNo = 1:length(EbNoVec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ofdm_symbols = 100;                       % Number of OFDM symbols to transmit
packet_length = n_ofdm_symbols*k*Ndata;   % Num de bits per packet

x = randint(1,packet_length);                   % Generate random data

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
% Transmit the signal (Raileigh, AWGN and add time/phase/freq noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we must pass the OFDMA signal on to a Raleigh channel
% Define the delay's vector for a CODIT Macro Model
% taken from "System-Level Simulation of a Third Generation WCDMA Wireless
% Geolocation Network" by Sanem Kabadayi.
%delays = [100e-9 200e-9 500e-9 600e-9 850e-9 900e-9 1050e-9 1350e-9 1450e-9 1500e-9];
%powers = [-3.2 -5.0 -4.5 -3.6 -3.9 -0.0 -3.0 -1.2 -5.0 -3.5];
delays = [Ts/50 Ts/20];
powers = [0 -3.2];
rayleigh = rayleighchan(1/Fs,40,delays,powers);
xrayleigh = filter(rayleigh,xrf);

% Add AWGN noise to the RF signal
ebno = EbNoVec(idxEbNo);
snr = ebno + 10*log10(k);
xnoisy = awgn(xrf, snr, 'measured');

yrf = xnoisy;
%yrf=xrf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S/P conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_sig_len = length(yrf);
cp_size = Nfft*G;
n_data_syms = floor(data_sig_len/(cp_size+Nfft));

% Cut to multiples of symbol period (Ts)
ycptime = reshape(yrf, (cp_size+Nfft), n_data_syms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the Cyclic Prefix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ycptime(1:cp_size,:) = [];
ytime = ycptime;

% Perform the FFT
yofdm = fft(ytime,Nfft);

% Extract data carriers
ymod = yofdm(UsedSubcIndex,:);
ymod = ymod(1:end);

% Extract pilot carriers
ypilot = yofdm(PilotSubcIndex,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Estimate (Known)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%channel_estimate = fft([zeros(size(rayleigh,1), abs(0)) rayleigh], Nfft, 2);
%channel_estimate = channel_estimate(:,UsedSubcIdx).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QPSK Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ysym = pskdemod(ymod,M);    

ybit = de2bi(ysym, 'left-msb');
y = reshape(ybit.',1,prod(size(ybit)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[number_of_errors(idxEbNo),bit_error_rate(idxEbNo)] = biterr(x,y);

end

% Plot the signal without and with AWGN noise.
figure(1);
plot(real(xrf(100:600)),'r');       % RF signal without AWGN
hold on;    
plot(real(yrf(100:600)),'b');       % RF signal with AWGN
hold off;
title('RF Signal Trasnmitted over AWGN channel');
xlabel('Amplitud'); ylabel('Samples');

% Plot power spectrum of OFDM signal
figure(2);
subplot(211);
[P,F] = pwelch(xrf,[],[],Nfft,Fs);
semilogy(F,P);
title('Tx OFDM Spectrum');
subplot(212);
[P,F] = pwelch(yrf,[],[],Nfft,Fs);
semilogy(F,P);
title('Tx OFDM Spectrum with AWGN');


% Plot BER performance
figure(3);
awgn_qpsk_theory = berawgn(EbNoVec, 'psk', 4, 'nondiff');

semilogy(EbNoVec, awgn_qpsk_theory, 'b-');
hold on;
semilogy(EbNoVec,bit_error_rate(:), 'r:*');

xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate');
grid on; drawnow;

title('Performance of 1024-FFT OFDMA-PHY');


