%
%   File:      ofdm_test1.m
%   Author(s): Horacio Sanson

%   Description:
%     Basic OFDMA-PHY tranceiver simulation.
%     Based on 128-FFT OFDMA-PHY (802.16e)
%     This similation uses Coded OFDM (COFDM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDMA-PHY Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Primitive Parameters (8.4.2.3)
Bw = 1.25e6;                                % Nominal Bandwidth 1.75MHz
n = 8/7;                                    % Oversampling Factor
G = 1/8;                                    % CP time to Tb factor
Nfft = 128;                                 % FFT Size

% Derived Parameters (8.4.2.4)
Fs = floor(n*Bw/8000)*8000;                 % Sampling Freq
deltaF = Fs/Nfft;                           % Subcarrier Spacing
Tb = 1/deltaF;                              % Useful symbol time
Tg = G*Tb;                                  % Guard Interval
Ts = Tb + Tg;                               % OFDMA Symbol time
T = Tb/Nfft;                                % Sampling Time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/2 Rate Convolutional Encoder Parameters  (8.4.9.2.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trellis_struct = poly2trellis( [7], [171,133]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Puncturing Parameters (Table 319)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code_rate_12 = 1/2;
punc_patt_12 = [1 1];                   % Puncturing pattern for 1/2 code rate
code_rate_23 = 2/3;         
punc_patt_23 = [1 0 1 1];               % Puncturing pattern for 2/3 code rate
code_rate_34 = 3/4;
punc_patt_34 = [1 0 1 1 1 0];           % Puncturing pattern for 3/4 code rate

code_rate = code_rate_12;
punc_patt = punc_patt_12;                % Change this value to choose the code rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CTC Interleaving Parameters (Does not follow the standard)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrows = 5;                      % Use 5 shift registers
slope = 3;                      % Delays are 0, 3, 6, 9 and 12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Modulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 4;                          % Alphabet Size (QPSK)
k = log2(M);                    % Number of bits per symbol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pilot/Data Subcarrier Allocation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation for 128-FFT FUSC (Table 309c)
Nleft = 11;                     % Left guard subcarriers 
Nrigth = 10;                    % Rigth guard subcarriers
Nused = 107;                    % Used subcarriers (Ndata+Ndc+Npilot)
Ndc = 1;                        % DC subcarrier, index 64
Ndata = 96;                     % Assign all data subcarriers (to simplify)
Npilot = 10;                    % Num pilot subcarriers

UsedSubcIndex = [12:74 75 76:118]';
PilotSubcPatt = [1 13 25 37 40 49 61 73 85 97]';
DataSubcPatt = [2:12 14:24 26:36 38:39 41:48 50:60 62:63 65:72 74:84 86:96 98:107]';
PilotSubcIndex = [12 24 36 48 51 60 72 84 96 108];
DataSubcIndex = [13:23 25:35 37:47 49:50 52:59 61:71 73:74 76:83 85:95 97:107 109:118]';

%%%%%%%%%%%%%%%%%%%%
% SIMULATION START %
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In OFDMA-PHY each user gets 48 subcarriers (1 subchannel) for transmit 
% in DL. We create packets with a size to fill 4 symbols of 48 subcarriers 
% each.

n_packets = 1;                                  % Num packets to transmit
n_ofdm_symbols = 4;                             % Number of OFDM symbols per packet
packet_length = n_ofdm_symbols*k*Ndata*code_rate;   % Num de bits per packet

x = randint(1,packet_length);                   % Generate random data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutional Encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcoded = convenc(x,trellis_struct);         % Conv Encode the data

% Puncture the coded data to change the code rate
pattern = repmat(punc_patt, 1, floor(size(xcoded,2)/size(punc_patt,2)));
rem_idx = find(pattern == 0);   % Get the indices that must be removed
xcoded(rem_idx) = [];           % Remove the indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CTC Interleaving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: No interleaving is implemented since it adds some delay to
% the signal. When I learn how to deal with the delay I will add it.
xinter = xcoded;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QPSK Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mapping = [0 1 3 2];        % Gray code mapping pattern for QPSK symbols
xbit = reshape(xinter,k,length(xinter)/k)';
xsym = bi2de(xbit, 'left-msb');
xmap = mapping(xsym+1);

xmod = pskmod(xmap,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Pilot Subcarriers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mod_syms = size(xmod,2);
n_ofdm_syms = n_mod_syms/Ndata;

% Create a new matrix and map xmod based on the DataSubcPatt
xofdm = zeros(Nused, n_ofdm_syms);
xofdm(DataSubcPatt,:) = reshape(xmod, Ndata, n_ofdm_syms);

% Fill the Pilot Subcarriers with the pattern (not implemented)
xofdm(PilotSubcPatt,:) = -1;   % Just fill them with -1's

% Serialize the data??
xpilot = xofdm(:).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to Time Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_symbols = size(xpilot,2)/Nused;

% Create a matrix to hold the subcarrier allocation
xtime = zeros(1, num_symbols*Nfft);

xpara = zeros(Nfft, num_symbols);
xpara(UsedSubcIndex,:) = reshape(xpilot, Nused, num_symbols);

xifft = ifft(xpara,Nfft);

xtime(1,:) = xifft(:).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the Cyclic Prefix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_symbols = size(xtime,2)/Nfft;
cp_size = Nfft*G;
tail_start = Nfft - cp_size + 1;
tail_end = Nfft;

xcptime = zeros(1,num_symbols*(cp_size+Nfft));
symbols = reshape(xtime(1,:), Nfft, num_symbols);
tmp_syms = [symbols(tail_start:tail_end,:); symbols];  % Add tail
%tmp_syms = [zeros(cp_size,num_symbols); symbols];  % Add zeros CP
xcptime(1,:) = tmp_syms(:).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit the signal (Raileigh, AWGN and add time/phase/freq noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xrf = xcptime;

yrf = xrf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return to Freq Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_sig_len = length(yrf);
cp_size = Nfft*G;
n_data_syms = floor(data_sig_len/(cp_size+Nfft));

% Cut to multiples of symbol period (Ts)
ytime = yrf(1:n_data_syms*(cp_size+Nfft));
yfft = reshape(ytime, (cp_size+Nfft), n_data_syms);

% Remove the guard interval
yfft(1:cp_size,:) = [];

% Perform the FFT
ypara = fft(yfft,Nfft);

% Extract data carriers
ymod = ypara(DataSubcIndex,:);
ymod = ymod(1:end);

% Extract pilot carriers
ypilot = ypara(PilotSubcIndex,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QPSK Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymap = pskdemod(ymod,M);    
[dummy demapping] = sort(mapping);
demapping = demapping -1;
ysym = demapping(ymap+1).';

ybit = de2bi(ysym, 'left-msb');
ycoded = reshape(ybit.',1,prod(size(ybit)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CTC Deinterleaving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yinter = ycoded;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Viterbi Decoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Depuncture the data by adding zeros
ydecoded = zeros(1,2*length(yinter)*code_rate);
pattern = repmat(punc_patt, 1, floor(size(ydecoded,2)/size(punc_patt,2)));
dat_idx = find(pattern == 1);
ydecoded(dat_idx)  = yinter;

% Decode using viterbi algorithm
y = vitdec(ydecoded,trellis_struct,96,'trunc','hard');




