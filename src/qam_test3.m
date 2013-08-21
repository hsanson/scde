%
%   File:      qam_test3.m
%   Author(s): Horacio Sanson
%
%   Description:
%     Emphyrical Simulation of QAM transmission/reception

M = 16;             % Alphabet Size
k = log2(M);        % Number of Bits per symbol
n = 3e4;            % Number of Bits to transmit
nsamp = 4;          % Oversampling rate

EbNo = [0:50];      % EbNo range in DB

% Initialize error vectors to store the simulated
% results
number_of_errors = zeros(1, length(EbNo));
bit_error_rate = zeros(1, length(EbNo));

% Initialize the root-raised cosine (Pulse Shaping Filter)
filtorder = 40; % Filter order
delay = filtorder/(nsamp*2); % Group delay (# of input samples)
rolloff = 0.5; % Rolloff factor of filter

% Define the delay's vector for a CODIT Macro Model
% taken from "System-Level Simulation of a Third Generation WCDMA Wireless
% Geolocation Network" by Sanem Kabadayi.
delays = [0 200e-9 500e-9];
powers = [0 -5.0 -4.5];
rayleigh = rayleighchan(1,0, delays, powers);
 
% Create a square root raised cosine filter.
rrcfilter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

% A. Define a vector for mapping bits to symbols using
% Gray coding. The vector is specific to the arrangement
% of points in a 16-QAM constellation.
mapping = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10].';

% Signal Source
x = randint(n,1);   % Random binary data
% B. Do ordinary binary-to-decimal mapping.
xsym = bi2de(reshape(x,k,length(x)/k).','left-msb');
% C. Map from binary coding to Gray coding.
xsym = mapping(xsym+1);

% Now transmit the signal throught a channel with
% different EbNo values.
for idxEbNo = 1:length(EbNo)
  
    % Modulate the signal
    ytx = qammod(xsym,M);  

    % Pulse Shape and oversample the signal
    ytx = rcosflt(ytx,1,nsamp,'filter',rrcfilter);
    
    % Pass the signal over a rayleigh fading channel
    ytx = filter(rayleigh,ytx);
   
    % Send the signal over an AWGN channel
    ebno = EbNo(idxEbNo);
    snr = ebno + 10*log10(k) - 10*log10(nsamp);
    ynoisy = awgn(ytx, snr, 'measured');
    
    % Received Signal and filter it with a square raised cosine filter
    yrx = rcosflt(ynoisy,1,nsamp,'Fs/filter',rrcfilter);
    yrx = downsample(yrx,nsamp); % Downsample.
    yrx = yrx(2*delay+1:end-2*delay); % Account for delay.
    
    % Demodulate the received signal
    ysym = qamdemod(yrx,M);
    
    % A. Define a vector that inverts the mapping operation.
    [dummy demapping] = sort(mapping);
    % Initially, demapping has values between 1 and M.
    % Subtract 1 to obtain values between 0 and M-1.
    demapping = demapping - 1;

    % B. Map between Gray and binary coding.
    ysym = demapping(ysym+1);

    % C. Do ordinary decimal-to-binary mapping.
    y = de2bi(ysym,'left-msb');
    % Convert z from a matrix to a vector.
    y = reshape(y.',prod(size(y)),1);
    
    % Calculate the bit error rate
    [number_of_errors(idxEbNo),bit_error_rate(idxEbNo)] = biterr(x,y);
end

figure(2);
ebno = 10.^(EbNo./10);  % Convert DB to ratio
fading_qam16_theory = (3/8).*(1-(1./sqrt(1+5./(2.*ebno))));

semilogy(EbNo, fading_qam16_theory, 'b-');
hold on;
semilogy(EbNo, bit_error_rate(:), 'b:*');

xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate');
grid on; drawnow;

title('Performance of QAM under Rayleigh Channel');


    
