
% DSP Handbook

% Example 17.4 Separation based on cycles
% Consider x(n) = s1(n)exp(j(w1n+theta1)) + s2(n)exp(j(w2n+theta2)) + v(n) 
% with s1(n), s2(n), and v(n) zero-mean Gaussian, and mutually independent.

% Let s1(n) be MA(3) with parameters [1, 0.2, 0.3, 0.5] and variance
% alpha1^2 = 1.38
% Let s2(n) be AR(1) with parameters [1, -0.5] and variance alpha2^2 = 2
% Let noise v(n) be MA(1) with parameters [1, 0.5] and variance alphav^2 =
% 1.25.

% Frequencies and phases are (w1, theta1) = (-0.5,0.6), (w2, theta2) =
% (1,1.8).


N = 1024;           % Number of samples
fftlength = 2048;   % FFT size
fs = 1;             % 1Hz
Ts = 1/fs;          % Sampling Period
n = [1:Ts:N];       % Time vector


% Create original signal
w1 = -0.5;                                  % s1 freq in radians
w2 = 1;                                     % s2 freq in radians
theta1 = 0.6;                               % s1 phase
theta2 = 1.8;                               % s2 phase

s1 = exp(i*(w1.*n + theta1));               % s1 signal
s2 = exp(i*(w2.*n + theta2));               % s2 signal

x = s1 + s2;                                % Original Signal

[Sx,fx] = pwelch(real(x),[],[],fftlength,fs);  % Signal PSD

figure(1);
subplot(411);
plot(n,x);
title('Original Signal');
subplot(412);
semilogy(fx*2*pi,Sx);
title('Original Signal PSD');

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
alphav = 1.25;
polyv = [1, 0.5];
v = sqrt(alphav)*randn(1,length(n));
v = filter(polyv,1,v);          % MA(1) process

% Generate the noisy signal 
xnoisy = s1noise.*s1 + s2noise.*s2 + v;          
[Sxnoisy,fxnoisy] = pwelch(real(xnoisy),[],[],fftlength,fs);


% Obtain the Cyclic Mean
Ak = [-4:0.5:4];
Cxx = zeros(size(Ak));

for I = 1:length(Ak)
    Cxx(I) = 0;
    for J = 1:N
        Cxx(I) = Cxx(I) + xnoisy(J)*xnoisy(J)*exp(-i*Ak(I)*J);
    end
    Cxx(I) = abs(Cxx(I)/N);
end


subplot(413);
plot(n,xnoisy);
subplot(414);
semilogy(fxnoisy*2*pi,Sxnoisy);

figure(2);
plot(Ak,Cxx);

%p1 = real(fftshift(fft(s1,2048)));
%freq = [-2048/2:2048/2-1].*1/2*2048;
%plot(freq,p1);

