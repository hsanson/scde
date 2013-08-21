%
%   File:      powerspectra2.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/6
%
%   Description:
%
%     Theorical Power Spectral Density of an OFDM signal using my own derived
%     equation.
%

clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDMA-PHY Parameters (Wimax Standard)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Primitive Parameters (8.4.2.3)
Bw = 1e3;                                  % Nominal Bandwidth
n = 8/7;                                   % Oversampling Factor
G = 1/8;                                   % CP time to Tb factor
Nfft = 128;                                % FFT Size

% Derived Parameters (8.4.2.4)
Fs = (n*Bw/8000)*8000;                 	    % Sampling Freq
deltaF = Fs/Nfft;                           % Subcarrier Spacing
Tb = 1/deltaF;                              % Useful symbol time
Tg = G*Tb;                                  % Guard Interval
Ts = Tb + Tg;                               % OFDMA Symbol time
T = Tb/Nfft;                                % Sampling Time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a simple BPSK signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn =sign(rand(1,Nfft)-0.5);	% Create BPSK signal
sn([1:20]) = 0;			% Zero the left guard bands 
sn([Nfft-20:Nfft]) = 0;         % Zero the rigth guard bands 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectra for one carrier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [-Fs:Fs];
k = [-Nfft/2:Nfft/2-1];
fk = k./Ts;

S = zeros(length(f),length(fk));

for i = 1:length(fk)

	C = (Ts/sqrt(T))*sn(i);
	W = sinc((f-fk(i))*Ts);
	P = exp((-j*2*pi*(1/2-(Tg/Ts))).*(f-fk(i)));

	S(:,i) = C .* W .* P;
end

figure(1);
subplot(211);
plot(f,abs(S).^2);
%axis([-Fs/2 Fs/2 0 1.5]);

subplot(212);
plot(f,abs(sum(S,2)));
%axis([-Fs/2 Fs/2 0 10e1]);
