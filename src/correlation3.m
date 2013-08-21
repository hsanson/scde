%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Cross Correlation with added noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
fs = 100;
Ts = 1/fs;
n = [0:Ts:150];
fftlength=2*length(n);
delay = 300;                % delay in samples (delay*Ts = sec)

w1 = -0.5;                                  % s1 freq in radians
w2 = 1;                                     % s2 freq in radians
theta1 = 0.6;                               % s1 phase
theta2 = 1.8;                               % s2 phase

s1 = exp(i*(w1.*n + theta1));               % s1 signal
s2 = exp(i*(w2.*n + theta2));               % s2 signal

% Multiplicative Noise (Rayleigh)
alpha1 = 1.38;
poly1 = [1, 0.2, 0.3, 0.5];
s1noise = sqrt(alpha1)*randn(1,length(n));
s1noise = filter(poly1,1,s1noise);        % MA(3) process

alpha2 = 2;
poly2 = [1, -0.5];
s2noise = sqrt(alpha2)*randn(1,length(n)); 
s2noise = filter(poly2,1,s2noise);        % AR(1) process

% Additive Noise (AWGN)
alphav1 = 1.28;
polyv1 = [1, 0.5];
v1 = sqrt(alphav1)*randn(1,length(n));
v1 = filter(polyv1,1,v1);          % MA(1) process

alphav2 = 0.31;
polyv2 = [1, 1.2];
v2 = sqrt(alphav2)*randn(1,length(n));
v2 = filter(polyv2,1,v2);          % MA(1) process

% Generate the noisy signal 
x = s1 + s2;
%xnoisy = s1noise.*x + v1;  
xnoisy = x + v1;

y = filter([1 1], 1, x);  % Create delayed version of x
y = [zeros(1,delay) y];
y = y(1:length(x));
%ynoisy = s2noise.*y + v2;
ynoisy = y + v2;

% Get power spectra of x, y and cross power spectra
[Sxx, Fxx] = pwelch(xnoisy,[],[],fftlength,fs);
[Syy, Fyy] = pwelch(ynoisy,[],[],fftlength,fs);
[Sxy, Fxy] = cpsd(xnoisy,ynoisy,[],[],fftlength,fs);

figure(1);
subplot(521); plot(n,x); title('Original signal x');
subplot(522); plot(n,xnoisy); title('Noisy signal x');
subplot(523); plot(n,y); title('Delayed signal y');
subplot(524); plot(n,ynoisy); title('Delayed Noisy signal y');
subplot(525); plot(Fxx,Sxx); title('Signal x PSD (Sxx)');
subplot(526); plot(Fyy,Syy); title('Delayed Signal y PSD (Syy)');
subplot(527); plot(Fxy,Sxy); title('Cross PSD Sxy');

% Obtain generalized cross correlation 
[corr1, t1] = GCC('unfiltered',Sxx, Syy, Sxy,fs,length(n));
[corr2, t2] = GCC('scot', Sxx, Syy, Sxy,fs,length(n));
[corr3, t3] = GCC('roth', Sxx, Syy, Sxy,fs,length(n));
[corr4, t4] = GCC('phat', Sxx, Syy, Sxy,fs,length(n));
[corr5, t5] = GCC('cps-m', Sxx, Syy, Sxy,fs,length(n));
[corr6, t6] = GCC('ht', Sxx, Syy, Sxy,fs,length(n));

figure(2);
subplot(321);
plot(t1,corr1);
axis([-5 5 -1 2]);
[m,k] = max(corr1);
k = t1(k);
title('Unfiltered GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));

subplot(322);
plot(t2,corr2);
axis([-5 5 -1 2]);
[m,k] = max(corr2);
k = t2(k);
title('SCOT GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));

subplot(323);
plot(t3,corr3);
axis([-5 5 -1 2]);
[m,k] = max(corr3);
k = t3(k);
title('ROTH GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));

subplot(324);
plot(t4,corr4);
axis([-5 5 -1 2]);
[m,k] = max(corr4);
k = t4(k);
title('PHAT GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));

subplot(325);
plot(t5,corr5);
axis([-5 5 -1 2]);
[m,k] = max(corr5);
k = t5(k);
title('cps-m GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));

subplot(326);
plot(t6,corr6);
axis([-5 5 -1 2]);
[m,k] = max(corr6);
k = t6(k);
title('HT GCC');
text(0.7,1.5,strcat(num2str(k),'sec'));


