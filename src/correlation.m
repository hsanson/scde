
% Proakis DSP

% Example 2.6.3 
% Shows that using correlation one can extract the period of an unknown
% signal even in the presence of noise. Looking at the correlation Ryy we
% see a period of N=10 samples even tought we contaminated the signal.

fc = 1/10;
fs = 5;
Ts = 1/fs;
n = [0:Ts:50];


xt = sin(2*pi*fc*n);

dB = 1;
Eb = sum(abs(xt).^2)/length(xt);
No = Eb / (10^(dB/10));

noise_var = sqrt(No/2);


nt = noise_var * randn(1, length(xt));

yt = xt + nt;

figure(1);
subplot(411);
stem(n,xt);
title('x[t] = sin 0.2\pi for 0<=n<=50');

subplot(412);
stem(n,nt);
title('nt[t] noise at 1dB SNR');

subplot(413);
stem(n,yt);
title('y[t] = x[t] + n[t]');

subplot(414);
len = length(yt);
stem([-len:len-2],xcorr(yt));
title('Ryy');
