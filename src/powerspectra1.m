%
%   File:      powerspectra1.m
%   Author(s): Horacio Sanson
%
%   Description:
%
%     Small Matlab program to show the different methods to estimate POWER
%     SPECTRUM DENSITIES.
%     This demo also demostrates the Parseval theorem. The average power of the signal
%     must be approximately equal to the average power of the spectrum estimates.
%     The matlab PSD must be scaled by a factor of Nfft/fs for the Parseval theorem to hold.
%
%   PSD Scaling:
%    http://www.mathworks.com/matlabcentral/newsreader/view_thread/149452#376001
%    http://groups.google.com/group/comp.soft-sys.matlab/browse_thread/thread/e8aaf649b67a9663/2222327db2ea7f51?lnk=st&q=greg-heath+power-spectrum-scaling&rnum=2#2222327db2ea7f51
%    http://www.dsprelated.com/showmessage/37758/1.php

fc = 25;                     % Signal freq
fs = 100;                    % Sampling freq
ts = 1/fs;                   % Sampling period
N = 100;                     % Number of samples
Nfft = 1024;                 % ffT size
scaling = Nfft/fs;           % Scaling factor needed to correct Matlab estimates
                             % for the Parseval theorem to hold.

n = [0:N-1]*ts;                     % Time vector
f = [-Nfft/2:Nfft/2-1]*fs/Nfft;     % freq vectors (-fs/2,fs/2) 

x = sin(2*pi*fc*n);                 % Generate signal
Px = sum(x.^2)/length(x)            % Obtain average power of signal

figure(1);

subplot(411);
Pxx = periodogram(x,[],Nfft,fs,'twosided');
Pxx = Pxx/scaling;
Pperio = sum(Pxx)   
plot(f,fftshift(Pxx));
title('PSD Estimate using Periodogram Method');

subplot(412);
Pxx = pwelch(x,[],[],Nfft,fs,'twosided');
Pxx = Pxx/scaling;
Pwelch = sum(Pxx)
plot(f,fftshift(Pxx));
title('PSD Estimate using Welch Method');

subplot(413);
Pxx = pmtm(x,[],Nfft,fs,'twosided');
Pxx = Pxx/scaling;
Pmtm = sum(Pxx)
plot(f,fftshift(Pxx));
title('PSD Estimate using  Thomson multitaper Method');

subplot(414);
Pxx = cpsd(x,x,[],[],Nfft,fs,'twosided');
Pxx = Pxx/scaling;
Pcpsd = sum(Pxx)
plot(f,fftshift(Pxx));
title('Cross PSD Estimate using Welch Method');
