%
%   File:      tde0.m
%   Author(s): Horacio Sanson
%
%   Description:
%     Cramer Rao Lower Bound on the variance of Time Delay Estimation of
%     Wideband signals. Taken from "An Overview on the Time Delay Estimate in
%     Active and Passive Systems for Target Localization" by Quazi.

% Initialize Matlab workspace
clear; clc;

SNR = [-50:50];		% SNR values in DB
snr = 10.^(SNR./10);	% SNR values in ordinary ratio

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

fo = 0;					            % Center frequency in Hz
f2 = fo + Bw/2;
f1 = fo - Bw/2;
T = Nfft*Ts;			            % Observation time in segs

% For low SNR (SNR << 1) in passive systems the CRLB on TDE is given by
var_low_snr = (3/(8*pi^2*T))^1/2 .* 1./snr .* 1/sqrt(f2^3-f1^3);

% For high SNR (SNR >> 1) in passive systems the CRLB on TDE is given by
var_high_snr = (3/(8*pi^2*T))^1/2 .* 1./sqrt(snr) .* 1/sqrt(f2^3-f1^3);

figure(1);

plot(SNR,var_low_snr);
hold on;
plot(SNR,var_high_snr, 'r');
legend('Low SNR', 'High SNR');


%% Numerical results %%

trials = 100;			% Number of trials to obtain variance
delay = 60;			% Delay in samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a simple BPSK signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn =sign(rand(1,Nfft)-0.5);	% Create BPSK signal
sn([1:90]) = 0;			% Zero the left guard bands (92 for 1024 Nfft)
sn([Nfft-01:Nfft]) = 0;         % Zero the rigth guard bands (91 for 1024 Nfft)

tde_variances = zeros(2, length(SNR));		% Store the variances for each GCC method

res = (Fs*Nfft)/(2*length(sn));

for idx = 1:length(snr)

	Eb = sum(abs(sn).^2)/length(sn);
	No = Eb/snr(idx);
	noise_var = sqrt(No./2);

	tde_trials = zeros(2, trials);

	for trial = 1:trials

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Create delayed signal
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		sn_delayed = filter([1 1], 1, sn);  % Create delayed version of sn
		sn_delayed = [zeros(1,delay) sn_delayed];
		sn_delayed = sn_delayed(1:length(sn));

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Add AWGN noise
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		n1 = noise_var * randn(1, length(sn));
		n2 = noise_var * randn(1, length(sn));

		sn_noisy = sn + n1;
		sn_delayed_noisy = sn_delayed + n2;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Estimate the spectra.
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[Sxx, Fxx] = pwelch(sn_noisy,[],[],Nfft, Fs, 'two-sided');
		[Syy, Fyy] = pwelch(sn_delayed_noisy,[],[],Nfft, Fs, 'two-sided');
		[Sxy, Fxy] = cpsd(sn_noisy,sn_delayed_noisy,[],[], Nfft, Fs, 'two-sided');

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Obtain the Generalized Correlation
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[corr1, t1] = GCC('unfiltered',Sxx, Syy, Sxy, Fs, length(sn_noisy), Nfft);
		[corr2, t2] = GCC('ht', Sxx, Syy, Sxy, Fs, length(sn_noisy), Nfft);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Obtain the sample number of the correlation max value
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[m1,k1] = max(corr1); tde_trials(1, trial) = t1(k1);
		[m2,k2] = max(corr2); tde_trials(2, trial) = t1(k2);

	end

	tde_variances(1, idx) = sum(((delay/res) - tde_trials(1,:)).^2)/length(tde_trials(1,:));
	tde_variances(2, idx) = sum(((delay/res) - tde_trials(2,:)).^2)/length(tde_trials(2,:));

end

hold off;
figure(2);
semilogy(SNR,tde_variances);
hold on;
semilogy(SNR,var_low_snr,'y+');
semilogy(SNR,var_high_snr,'c*');
title('Variance of Time Delay Estimates using GCC methods');
legend('Direct Cross Correlation', 'HT GCC', 'Low SNR CRLB', 'High SNR CRLB');

