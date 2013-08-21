% Time Delay estimation for OFDM signal.
% Uses a basic 1014-FFT OFDMA-PHY tranceiver
% Here we compare two different methods to obtain the correlation functions.
% One using the serialized rf signal and the other using averaging over sub carriers.

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

n_ofdm_symbols = 2;                        % Number of OFDM symbols to transmit
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

EbNo = [-15:10];      	% EbNo range in DB
trials = 10;		% Number of simulation runs per SNR value to obtain variance estimates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tde_variances = zeros(4, length(EbNo));		% Store the variances for each GCC method
tde_crlb = zeros(1,length(EbNo));
res = (Fs*Nfft)/(2*length(yrf));

% For each SNR value
for idxEbNo = 1:length(EbNo)

	tde_trials = zeros(4, trials);

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
		% S/P conversion
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		data_sig_len = length(xrf);
		cp_size = Nfft*G;
		n_data_syms = floor(data_sig_len/(cp_size+Nfft));

		% Cut to multiples of symbol period (Ts)
		xcptime = reshape(xnoisy, (cp_size+Nfft), n_data_syms);

		% Cut to multiples of symbol period (Ts)
		ycptime = reshape(ynoisy, (cp_size+Nfft), n_data_syms);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Estimate the spectra.
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% For the parallel case
		Sxxp = zeros(Nfft, size(xcptime,2));
		Syyp = zeros(Nfft, size(xcptime,2));
		Sxyp = zeros(Nfft, size(xcptime,2));

		% For each subcarrier
		for carrier=1:size(ycptime,2)
			Sxxp(:,carrier) = pwelch(xcptime(:,carrier),[],[],Nfft, 1/Ts);
			Syyp(:,carrier) = pwelch(ycptime(:,carrier),[],[],Nfft, 1/Ts);
			Sxyp(:,carrier) = cpsd(xcptime(:,carrier),ycptime(:,carrier),[],[], Nfft, 1/Ts);
		end

		% For the serial case
		Sxxs = pwelch(xnoisy,[],[],Nfft, Fs);
		Syys = pwelch(ynoisy,[],[],Nfft, Fs);
		Sxys = cpsd(xnoisy,ynoisy,[],[], Nfft, Fs);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Obtain the Generalized Correlation
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[corr1p, t1] = GCC('unfiltered',Sxxp, Syyp, Sxyp, Fs, length(xnoisy), Nfft);
		[corr2p, t2] = GCC('ht', Sxxp, Syyp, Sxyp, Fs, length(xnoisy), Nfft);
		[corr1s, t3] = GCC('unfiltered', Sxxs, Syys, Sxys, Fs, length(xnoisy), Nfft);
		[corr2s, t4] = GCC('ht', Sxxs, Syys, Sxys, Fs, length(xnoisy), Nfft);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Average the correlation values (only for parallel case)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		corr1p = mean(corr1p,2);
		corr2p = mean(corr2p,2);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Obtain the sample number of the correlation max value
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[m1,k1] = max(corr1p); tde_trials(1, trial) = t1(k1);
		[m2,k2] = max(corr2p); tde_trials(2, trial) = t2(k2);
		[m3,k3] = max(corr1s); tde_trials(3, trial) = t3(k3);
		[m4,k4] = max(corr2s); tde_trials(4, trial) = t4(k4);

	end

	tde_variances(1, idxEbNo) = sum(((delay/res) - tde_trials(1,:)).^2)/length(tde_trials(1,:));
	tde_variances(2, idxEbNo) = sum(((delay/res) - tde_trials(2,:)).^2)/length(tde_trials(2,:));
	tde_variances(3, idxEbNo) = sum(((delay/res) - tde_trials(3,:)).^2)/length(tde_trials(3,:));
	tde_variances(4, idxEbNo) = sum(((delay/res) - tde_trials(4,:)).^2)/length(tde_trials(4,:));

	%a = (3/(8*pi^2*length(x)))^1/2
	%b = (1/sqrt(Bw^3 - 0^3))
	%c = (1/(10^(snr/10)))
	%crlb = a * c * b
	%tde_variances(5,idxEbNo) = (3/(8*pi^2*length(x)))^1/2  * (1/(10^(snr/10))) * (1/sqrt(Bw^3 - 0^3))
	%tde_crlb(1,idxEbNo) = crlb;
end

figure(1);
subplot(2,1,1);
semilogy(EbNo,tde_variances([1 3], :));
title('Variance of Time Delay Estimates Direct Correlation');
%hold on;
%plot(EbNo, tde_crlb(1,:));
legend('averaging', 'no averaging');

subplot(2,1,2);
semilogy(EbNo,tde_variances([2 4], :));
title('Hannan Thompson Generalized Cross Correlation');
legend('averaging', 'no averaging');



