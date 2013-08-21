%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      bpsk_test1.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%     Emphyrical Simulation of BPSK transmission/reception
clear; clc;

M = 2;              % Alphabet Size
k = log2(M);        % Number of Bits per symbol
n = 3e4;            % Number of Bits to transmit
nsamp = 2;          % Oversampling rate

EbNo = [-10:5];      % EbNo range to simulate

% Initialize error vectors to store the simulated
% results
number_of_errors = zeros(1, length(EbNo));
bit_error_rate = zeros(1, length(EbNo));

% Initialize the root-raised cosine (Pulse Shaping Filter)
filtorder = 40; % Filter order
delay = filtorder/(nsamp*2); % Group delay (# of input samples)
rolloff = 0.5; % Rolloff factor of filter

% Create a square root raised cosine filter.
rrcfilter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

% Signal Source
x = randint(n,1);   % Random binary data

% Now transmit the signal throught a channel with
% different EbNo values.
for idxEbNo = 1:length(EbNo)

    % Modulate the signal
    ytx = pskmod(x,M);

    % Pulse Shape and oversample the signal
    ytx = rcosflt(ytx,1,nsamp,'filter',rrcfilter);

    % Send the signal over an AWGN channel
    ebno = EbNo(idxEbNo);
    snr = ebno + 10*log10(k) - 10*log10(nsamp);
    ynoisy = awgn(ytx, snr, 'measured');

    % Received Signal and filter it with a square raised cosine filter
    yrx = rcosflt(ynoisy,1,nsamp,'Fs/filter',rrcfilter);
    yrx = downsample(yrx,nsamp); % Downsample.
    yrx = yrx(2*delay+1:end-2*delay); % Account for delay.

    % Demodulate the received signal
    y = pskdemod(yrx,M);

    % Calculate the bit error rate
    [number_of_errors(idxEbNo),bit_error_rate(idxEbNo)] = biterr(x,y);
end

figure(1);
awgn_bpsk_theory = berawgn(EbNo, 'psk', 2, 'nondiff');

semilogy(EbNo, awgn_bpsk_theory, 'b-');
hold on;
semilogy(EbNo, bit_error_rate(:), 'r:*');

xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate');
grid on; drawnow;

title('Performance of BPSK');
legend('Theorical', 'Simulated');

