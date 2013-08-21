%
%   File:      powerspectra9.m
%   Author(s): Horacio Sanson
%   Revision : 2007/10/24
%
%   Description:
%     SCD estimation example taken from the following forum thread:
%     http://bbs.bc-cn.net/TopicOther.asp?t=5&BoardID=216&id=76056

Rb=2;     % Symbol Rate
Ts=1/Rb;  % Symbole time
L=2^6;    % Oversampling rate for pulse shaping (64)
num=16;   % Number of BPSK symbols
N=num*L;  % Total number of samples (after pulse shaping)
f0=16;    % Carrier Frequency
dt=Ts/L;  %抽样时间间隔
df=1/(dt*N);   %频率间隔

Bs=df*N/2;   %带宽
fs=Bs*2;
f=df/2-Bs:df:Bs;   %抽样频率
t=dt/2:dt:num*Ts;  %抽样时间

% random BPSK signal generation
a=sign(randn(1,num));

% Pulse shaping, this looks to be the same as rectpulse(a,L);
for rr=1:num
    I((rr-1)*L+1:rr*L)=a(rr);
end

% Carrier Modulate
s=I.*cos(2*pi*f0*t);


M=2;
Fs=fs;
N = (M*Fs)/df;
N = pow2 (nextpow2(N)); % windowing record for FFT
%N=; %设置FFT点数为1024
X = fft(s,N);   % fft of the truncated (or zero padded) time series
X = fftshift(X);% shift components of fft
Xc = conj(X);                 % precompute the complex conjugate vector
S = zeros (N,N);              % size of the Spectral Correlation Density matrix
f = zeros (N,N);              % size of the frequency matrix;
alfa = zeros (N,N);           % size of the cycle frequency matrix
F = Fs/(2*N);                 % precompute constants -  F = Fs/(2*N);
G = Fs/N;                     % precompute constants -  G = Fs/N;
m = -M/2+1:M/2;               % set frequency smoothing window index
for k = 1:N                              % fix k
    % computes vectors of f and alfa,
    % store frequency and cycle frequency data for given k.

    k1 = 1:N;
    f(k,k1) = F*(k+k1-1) - Fs/2;          % Computes f values and shift them to center in zero (f = (K+L)/2N) [1]
    alfa(k,k1) = G*(k-k1 + N-1) - Fs;       % Computes alfa values and shift them to center in zero (alfa = (K-L)/N) [1]

    for k1 = 1:N %fix k1 = J
        %calculate X(K+m) & conj (X(J+m)) for arguments of X(1:N) only
        B = max(1-k, 1-k1);          % Largest min of 1 <= (K+m)| (J+m) <= N
        A = min (N-k, N-k1);         % Smallest max of 1 <= (K+m)| (J+m) <= N
        n = m((m<=A) & (m>=B));   %fix the index out of range problem by
                                                   % truncating the window
        if isempty(n)
            S(k,k1) = 0;
        else
            p = k+n; q = k1+n;
            Y = X(p).*Xc(q);
            S(k,k1) = sum(Y);
        end
    end
end
S = abs(S./max(max(S)));% normalize output matrix
% figure(1);
% contour (alfa, f, S); grid;
figure(1);
mesh(alfa, f, S);
ylabel('Frequency (Hz)');
xlabel('Cycle frequency (Hz)');
title('BPSK循环谱三维图');
