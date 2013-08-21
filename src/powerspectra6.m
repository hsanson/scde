%
%   File:      powerspectra5.m
%   Author(s): Horacio Sanson

%   Description:
%
%     Power spectral density estimation example using different
%     covariance estimators. The methods for covariance estimation
%     we use here are those supported by the corrmtx Matlab function.
%     See corr_matlab.m for details.
%
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n         = 512;                       % Number of observed samples
A = [ 1 1 1 1 1 ];                     % Amplitudes of each sinusoid
f = [0.05;0.06;0.07;0.08;0.2];         % Frequencies of each sinusoid
fs        = 1;                         % Sampling frequency
variance  = 0.3;                       % Noise variance of the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nT        = [1:1/fs:n];             % Time axis
fftlength = 2^nextpow2(n);          % FFT length

% Cool trick that generates a signal composed of the sum of
% several sinusoids + white noise. (Matlab Specral Analysis Help)
x = A*cos(2*pi*f*nT) + variance * randn(size(nT));

M = length(x)/4;

[S,k] = capon(x,x,M,'autocorrelation');
S = 10*log10(S/max(S));

figure(1);
plot(k*fs,S, 'LineWidth', 1, 'Color', [0.6 0.7 1]);
title('Estimated PSD via CAPON method');
xlabel('freq (Hz)');
ylabel('PSD (dB)');

hold on;
[S,k] = capon(x,x,M,'prewindowed');
S = 10*log10(S/max(S));
plot(k*fs,S, 'LineWidth', 1, 'Color', [0.5 1.0 0.5]);

[S,k] = capon(x,x,M,'postwindowed');
S = 10*log10(S/max(S));
plot(k*fs,S, 'LineWidth', 1, 'Color', [1.0 1.0 0.3]);

[S,k] = capon(x,x,M,'covariance');
S = 10*log10(S/max(S));
plot(k*fs,S, 'LineWidth', 1, 'Color', [1.0 0.0 1.0]);

[S,k] = capon(x,x,M,'modified');
S = 10*log10(S/max(S));
plot(k*fs,S, 'LineWidth', 1, 'Color', [0.0 1.0 1.0]);
legend('autocorrelation','prewindowed', 'postwindowed', 'covariance', 'modified (Forward-Backward)');
