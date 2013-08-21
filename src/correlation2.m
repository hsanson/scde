%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Cross Correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 100;
Ts = 1/fs;
n = [0:Ts:20];
fftlength=2*length(n);
delay = 51;                % delay in samples (delay*Ts = sec)

x = randn(1,length(n));
%x = sin(2*pi*1*n);

y = filter([1 1], 1, x);  % Create delayed version of x
y = [zeros(1,delay) y];
y = y(1:length(x));

% Get power spectra of x, y and cross power spectra
[Sxx, Fxx] = pwelch(x,[],[],fftlength,fs);
[Syy, Fyy] = pwelch(y,[],[],fftlength,fs);
[Sxy, Fxy] = cpsd(x,y,[],[],fftlength,fs);

%figure(1);
%subplot(511);
%plot(n,x);
%subplot(512);
%plot(Fxx,Sxx);
%subplot(513);
%plot(n,y);
%subplot(514);
%plot(Fyy,Syy);
%subplot(515);
%plot(Fxy,Sxy);

Ts*delay

% Obtain cross correlation 
[corr1, t1] = GCC('unfiltered',Sxx, Syy, Sxy,fs,length(n));
[corr2, t2] = GCC('scot', Sxx, Syy, Sxy,fs,length(n));
[corr3, t3] = GCC('roth', Sxx, Syy, Sxy,fs,length(n));
[corr4, t4] = GCC('phat', Sxx, Syy, Sxy,fs,length(n));
[corr5, t5] = GCC('cps-m', Sxx, Syy, Sxy,fs,length(n));
[corr6, t6] = GCC('ht', Sxx, Syy, Sxy,fs,length(n));

figure(2);
subplot(321);
plot(t1,corr1);
axis([-25 25 -1 1]);
[m,k] = max(corr1);
k = t1(k);
title('Unfiltered GCC');
text(10,0.5,strcat(num2str(k),'sec'));

subplot(322);
plot(t2,corr2);
axis([-25 25 -1 1]);
[m,k] = max(corr2);
k = t2(k);
title('SCOT GCC');
text(10,0.5,strcat(num2str(k),'sec'));

subplot(323);
plot(t3,corr3);
axis([-25 25 -1 1]);
[m,k] = max(corr3);
k = t3(k);
title('ROTH GCC');
text(10,0.5,strcat(num2str(k),'sec'));

subplot(324);
plot(t4,corr4);
axis([-25 25 -1 1]);
[m,k] = max(corr4);
k = t4(k);
title('PHAT GCC');
text(10,0.5,strcat(num2str(k),'sec'));

subplot(325);
plot(t5,corr5);
axis([-25 25 -1 1]);
[m,k] = max(corr5);
k = t5(k);
title('cps-m GCC');
text(10,0.5,strcat(num2str(k),'sec'));

subplot(326);
plot(t6,corr6);
axis([-25 25 -1 1]);
[m,k] = max(corr6);
k = t6(k);
title('HT GCC');
text(10,0.5,strcat(num2str(k),'sec'));


