function [G,t,R] = GCC(m,Pxx,Pyy,Pxy,Fs,frame,N)
%GCC
% Generalized Cross-Correlation with specified pre-whitening filter.
%
% Description: the GCC is computed with a pre-whitening filter onto the
%              cross-power spectrum in order to weight the magnitude value
%              against its SNR. The weighted CPS is used to obtain the
%              cross-correlation in the time domain with an inverse Fourier
%              transformation. The result is normalized between [-1 1].
%
% Usage      : G = GCC(m,Pxx,Pyy,Pxy)
%              [G,t] = GCC(m,Pxx,Pyy,Pxy,Fs,frame)
%              [G,t] = GCC(m,Pxx,Pyy,Pxy,Fs,frame,N)
%              [G,t,R] = GCC(...)
%
% m            pre-whitening filter type
%              ('Roth','SCOT','PHAT','CPS-m','HT','unfiltered')
% Pxx,Pyy,Pxy  power spectra of the input signal. If Pxx,Pyy,Pxy are
%              matrices each column represent the power spectrum for the
%              relative channel involved in the Pxy calculus:
%              (i.e. Pxx contains power spectrum for channel (1) (1) (2),
%              Pyy contains power spectrum for channel (2) (3) (3), so Pxy
%              contains the cross power spetrum for the couples (1,2) (1,3)
%              (2,3))
% Fs           sampling frequency (Hz)
% frame        length in samples of the input sequence from which were
%              extracted the power spectra. It's used to compute the time
%              ticks
% N            [optional] FFT length. By default Pxy length
%
% G            vector of GCC values
%              matrix of GCC values on columns when Pxy is a matrix
% t            [optional] vector with time ticks for the GCC
%              [optional] vector of frequencies GCC values
% R            matrix of frequencies GCC values where columns are GCC for
%              the respective channel
%
%
% Author     : Davide Renzi, d.renzi@infocom.uniroma1.it
%              INFOCOM Dept. University of Rome "La Sapienza"
% Version    : 2.0
% Date       : Rome, June 2005
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.


% check input arguments
if nargout==1 & nargin<4
    error('Wrong number of input arguments');
elseif nargout==2 & nargin<6
    error('Wrong number of input arguments to compute time ticks');
elseif nargin<7
    N = size(Pxy,1);
end


% number of channels
Nchn = size(Pxy,2);

% define pre-whitening filter
switch lower(m)
    case 'unfiltered'
        % Unfiltered Cross-Correlation (UCC)
        % -----------------------------------------------------------------
        % this processor doesn't filter the cross-power applying a sort of
        % time cross-correlation.
        W = ones(N,Nchn);

    case 'roth'
        % Roth filter
        % -----------------------------------------------------------------
        % this processor suppress frequency regions where the noise is
        % large than signals.
        if size(Pxx)~=size(Pxy)
            error('Roth filter: power spectra size must be the same');
        end
        W = ones(N,Nchn);
        W_Den = Pyy;
        for k=1:Nchn
            nonzero = find(W_Den(:,k));
            W(nonzero,k) = 1 ./ W_Den(nonzero,k);
        end

    case 'scot'
        % Smoothed Coherence Transform (SCOT)
        % -----------------------------------------------------------------
        % this processor exhibits the same spreading as the Roth processor.
        if size(Pxx)~=size(Pyy) | size(Pxx)~=size(Pxy) | size(Pyy)~=size(Pxy)
            error('SCOT filter: power spectra size must be the same');
        end
        W = ones(N,Nchn);
        W_Den = (Pxx .* Pyy) .^ 0.5;
        for k=1:Nchn
            nonzero = find(W_Den(:,k));
            W(nonzero,k) = 1 ./ W_Den(nonzero,k);
        end

    case 'phat'
        % Phase Transform (PHAT)
        % -----------------------------------------------------------------
        % ad hoc technique devoleped to assign a specified weight according
        % to the SNR.
        W = ones(N,Nchn);
        W_Den = abs(Pxy);
        for k=1:Nchn
            nonzero = find(W_Den(:,k));
            W(nonzero,k) = 1 ./ W_Den(nonzero,k);
        end

    case 'cps-m'
        % SCOT filter modified
        % -----------------------------------------------------------------
        % this processor computes the Cross Power Spectrum Density and
        % apply the SCOT filter with a power at the denominator to avoid
        % ambient reverberations that causes false peak detection.
        if size(Pxx)~=size(Pyy) | size(Pxx)~=size(Pxy) | size(Pyy)~=size(Pxy)
            error('CPS-modified: power spectra size must be the same');
        end
        W = ones(N,Nchn);
        factor = .75; % common value between .5 and 1
        W_Den = (Pxx .* Pyy) .^ factor;
        for k=1:Nchn
            nonzero = find(W_Den(:,k));
            W(nonzero,k) = 1 ./ W_Den(nonzero,k);
        end

    case 'ht'
        % Hannah and Thomson filter (HT)
        % -----------------------------------------------------------------
        % HT processor computes a PHAT transform weighting the phase
        % according to the strength of the coherence.
        if size(Pxx)~=size(Pyy) | size(Pxx)~=size(Pxy) | size(Pyy)~=size(Pxy)
            error('HT filter: power spectra size must be the same');
        end
        W = ones(N,Nchn);
        gamma = ones(N,Nchn);
        % coherence function evaluated along the frame
        gamma_Den = (Pxx .* Pyy) .^ 0.5;
        for k=1:Nchn
            nonzero = find(gamma_Den(:,k));
            gamma(nonzero,k) = abs(Pxy(nonzero,k) ./ gamma_Den(nonzero,k)).^2;
        end
        % HT filter
        W_Den = abs(Pxy) .* (1-gamma);
        for k=1:Nchn
            nonzero = find(W_Den(:,k));
            W(nonzero,k) = gamma(nonzero,k) ./ W_Den(nonzero,k);
        end

    otherwise error('Method not defined...');

end


% apply the filter
R = Pxy .* W;

% estimate the generalized cross-correlation (GCC)
G = fftshift(real(ifft(R)),1);
% NB: the real part is extracted to avoid the small undesidered imag. part
%     that sometimes compare during the inverse Fourier transform.

% normalize GCC
for k=1:Nchn
    G(:,k) = G(:,k)/max(abs(G(:,k)));
end

if nargout>1
    % calculate thick along time axis
    resolution = (Fs*N)/(2*frame);
    if mod(N,2)==0
        t = [-N/2 : (N/2)-1] / resolution;
    else
        t = [-(N-1)/2 : (N-1)/2] / resolution;
    end
end
