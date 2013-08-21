%
%   File:      wosa_test02.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/6

%   Description:
%       Example of SCD estimation of Sinusoid signal using WOSA for diferent filter sizes.
%   Notes:
%        - Tested with Matlab 2007a

clear; clc;
warning off all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters of the simulation. All changes
% must be done here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = 64;              % Number of observed samples
A       = 1;               % SNR of each sinusoid
fc      = 1.5;             % Frequencies of each sinusoid
fs      = 10;              % Sampling frequency
trials  = 10;              % Number of monte carlo trials
snr     = 0;              % SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters that are obtained based on the
% General Parameters above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts      = 1/fs;               % Sampling period
nT      = [0:N-1]*ts;         % Time axis
M       = [N/16 N/8 N/4 N/2]; % Filter lengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the sinusoidal signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = A*sin(2*pi*fc*nT);
y = awgn(y,snr,'measured');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the theorical PSD of the generated signal
% only works for one sinusoid in noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h   = zeros(N,N);           % Array to hold theorical value
f  = [-N/2:N/2-1]*fs/N;     % frequency axis, same for cyclic freq
fa  = [-N/2:N/2-1]*fs/N;     % frequency axis, same for cyclic freq
f1 = [-fliplr(fc); fc];     % frequency bins  
f2 = [-fliplr(2*fc); 2*fc]; % alpha bins (cyclic freq at -2fc and 2fc)
a1 = [fliplr(A) A];         % bins amplitude

for i = 1:length(f1)
    [m,ix] = min(abs(f-f1(i)));   % find index of the closest value to the frequency
    h(N/2,ix) = (a1(i)^2)/2;      % theorical PSD
end

for i = 1:length(f2)
    [m,ix] = min(abs(f-fa(i)));   % find index of the closest value to the frequency
    h(ix,N/2) = (a1(i)^2)/2;           % theorical PSD
end
h = h';
hdb = 10*log10(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the SCDs via different methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,f,a] = scd(y,N/16,'WOSA',fs);
figure(1);
subplot(221);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(221);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/16)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/8,'WOSA',fs);
figure(1);
subplot(222);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(222);
mesh(a,f,abs(S));
title('BPSK Spectral Correlation (WOSA N/8)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/4,'WOSA',fs);
figure(1);
subplot(223);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(223);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/4)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

[S,f,a] = scd(y,N/2,'WOSA',fs);
figure(1);
subplot(224);
contour(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);
figure(4);
subplot(224);
mesh(a,f,abs(S));
title('Sinusoid Spectral Correlation (WOSA N/2)','FontSize', 16);
xlabel('\alpha','FontSize', 16);
ylabel('freq (Hz)','FontSize', 16);
zlabel('PSD','FontSize', 16);

