%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      sampling.m
%   Author(s): Horacio Sanson
%
%   Description:
%     Simple example of a sin wave with spectrum

fs = 100;               % Sampling frequency in Hz -> number of samples per 1 seg
Ts = 1/fs;              % Sampling period
n = [1:1024]*Ts;          % Number of samples to process (5 sec)
fftlength = 1024;       % Size of the fft (the larger the better)


% Generate 2Hz sin signal
fc = 2;                % 2 Hz
x = sin(2*pi*fc.*n);   % Sinusoid signal

% Generate 3Hz sin signal
fc = 20;               % 10 Hz
y = sin(2*pi*fc.*n);   % Sinusoid signal

% Plot signal
figure(1);
subplot(311);
plot(n,x);
hold on;
plot(n,y,'r-');
title('Time Domain Signals');
hold off;

% Plot spectrum
% Obtain the x scale in Hz. Remember that FFT corresponds to th sampling
% rate fs for continuous frequencies and corresponds to 2*pi for discrete
% frequencies.
k = [-fftlength/2:fftlength/2-1].*fs/fftlength;


subplot(312);
plot(k,real(fftshift(fft(x,fftlength))));
hold on;
plot(k,real(fftshift(fft(y,fftlength))), 'r');
hold off;


p = abs(fft(x,fftlength)/length(x)/2).^2;
subplot(313);
plot(p(1:length(x)/2));
