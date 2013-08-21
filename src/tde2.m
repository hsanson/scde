%
%   File:      tde2.m
%   Author(s): Horacio Sanson
%
%   Description:
%     Time Delay estimation for OFDM signal.
%     Uses a basic 1014-FFT OFDMA-PHY tranceiver

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

%%%%%%%%%%%%%%%%%%%%
% SIMULATION START %
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ofdm_symbols = 1;                 % Number of OFDM symbols to transmit
nbits = n_ofdm_symbols*k*Ndata;     % Num de bits to transmit. We make sure they
				    % fit exactly on the number of OFDM symbols

x = randint(1,nbits);               % Generate random data

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
% xofdm(PilotSubcIndex,:) = -1;   % Just fill them with -1's

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
delay = 60;		      % delay in samples

yrf = filter([1 1], 1, xrf);  % Create delayed version of xrf
yrf = [zeros(1,delay) yrf];
yrf = yrf(1:length(xrf));

EbNo = [-50:10];      	% EbNo range in DB
trials = 100;		% Number of simulation runs per SNR value to obtain variance estimates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tde_variances = zeros(6, length(EbNo));		% Store the variances for each GCC method 
tde_mean_error = zeros(6, length(EbNo));	% Store the mean error of TDE

res = (Fs*Nfft)/(2*length(yrf));

% For each SNR value 
for idxEbNo = 1:length(EbNo)

	tde_trials = zeros(6, trials);

	EbNo(idxEbNo)				% Print current processing EbNo 

	for trial = 1:trials

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Transmit both signals throught an AWGN Channel 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		ebno = EbNo(idxEbNo);
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

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Obtain the sample number of the correlation max value
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[m1,k1] = max(corr1); tde_trials(1, trial) = t1(k1);
		[m2,k2] = max(corr2); tde_trials(2, trial) = t2(k2);
		[m3,k3] = max(corr3); tde_trials(3, trial) = t3(k3);
		[m4,k4] = max(corr4); tde_trials(4, trial) = t4(k4);
		[m5,k5] = max(corr5); tde_trials(5, trial) = t5(k5);
		[m6,k6] = max(corr6); tde_trials(6, trial) = t6(k6);

	end

	tde_variances(1, idxEbNo) = sum(((delay/res) - tde_trials(1,:)).^2)/length(tde_trials(1,:));
	tde_variances(2, idxEbNo) = sum(((delay/res) - tde_trials(2,:)).^2)/length(tde_trials(2,:));
	tde_variances(3, idxEbNo) = sum(((delay/res) - tde_trials(3,:)).^2)/length(tde_trials(3,:));
	tde_variances(4, idxEbNo) = sum(((delay/res) - tde_trials(4,:)).^2)/length(tde_trials(4,:));
	tde_variances(5, idxEbNo) = sum(((delay/res) - tde_trials(5,:)).^2)/length(tde_trials(5,:));
	tde_variances(6, idxEbNo) = sum(((delay/res) - tde_trials(6,:)).^2)/length(tde_trials(6,:));

end	

figure(1);
plot(EbNo,tde_variances);
title('Variance of Time Delay Estimates using GCC methods');
legend('Direct Cross Correlation', 'SCOT GCC', 'Roth GCC', 'Phat GCC', 'CPS-M GCC', 'HT GCC');

%figure(2);
%plot(EbNo,tde_mean_error);
%title('Mean Error of Time Delay Estimates using GCC methods');
%legend('Direct Cross Correlation', 'SCOT GCC', 'Roth GCC', 'Phat GCC', 'CPS-M GCC', 'HT GCC');
	
figure(3);
plot(EbNo,tde_variances([1 6],:));
title('Variance of Time Delay Estimates using GCC methods');
legend('Direct Cross Correlation', 'HT GCC');

%figure(4);
%plot(EbNo,tde_mean_error([1 6],:));
%title('Mean Error of Time Delay Estimates using GCC methods');
%legend('Direct Cross Correlation', 'HT GCC');
