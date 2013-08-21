%
%   File:      powerspectra11.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/17
%
%   Description:
%
%     Simple test to see the periodic nature of the DFT.
%     I generate a sin wave and plot it's FFT, then shift the signal
%     by multiplying it by an exponential and plot the FFT again. If we keep shifting
%     the sinusoidal component will wrap around the graph

clear; clc;

%% Create a Sinusoid signal

fc = 2;                        % Signal freq
fs = 10;                       % Sampling Freq
Ts = 1/fs;
N = 512;
n = [0:N-1]'*Ts;
variance = 0.1;
Nfft = 512;
shift = 5;                      % freq shift in Hz
sn_noisy = sin(2*pi*fc.*n) + rand(length(n),1)*variance;

%%Shift the signal


f = [-Nfft/2:Nfft/2-1]*fs/Nfft;

figure(1);
plot(f,fftshift(abs(real(fft(sn_noisy,Nfft)))));
hold on;

ex = exp(j*pi*shift*n);

plot(f,fftshift(abs(real(fft(sn_noisy.*ex,Nfft)))),'r^-.');

ex = exp(-j*pi*2*shift*n);

plot(f,fftshift(abs(real(fft(sn_noisy.*ex,Nfft)))),'g+:');

legend('Original', 'Shift 5Hz', 'Shift 10Hz');
