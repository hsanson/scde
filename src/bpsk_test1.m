%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      bpsk_test1.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%    Simple comparison between BPSK and QPSK modulation.

clear; clc;

EbNoVec = [0:10];

% Matlab accepts EbNo valies in DB but in most books the EbNo is shown as a
% ratio (http://www.transcrypt.com/download?id=7550).

ebno = 10.^(EbNoVec./10);   % Convert EbNo from DB to ratio

ber_awgn_bpsk = 1/2 .* erfc(sqrt(ebno));  % Obtain BER value

% Theorical BER as obtained with the berawgn Matlab function

% BER on AWGN for BPSK
ber_awgn_bpsk1 = berawgn(EbNoVec, 'psk', 2, 'nondiff');
% BER on AWGN for BPSK with differential phase
ber_awgn_bpsk2 = berawgn(EbNoVec, 'psk', 2, 'diff');
% BER on AWGN for QPSK
ber_awgn_qpsk1 = berawgn(EbNoVec, 'psk', 4, 'nondiff');
% BER on AWGN for QPSK with differential phase
ber_awgn_qpsk2 = berawgn(EbNoVec, 'psk', 4, 'diff');

figure(1);
semilogy(EbNoVec,ber_awgn_bpsk1, 'm');
hold on;
semilogy(EbNoVec,ber_awgn_bpsk2, 'r');
semilogy(EbNoVec,ber_awgn_qpsk1, 'g');
semilogy(EbNoVec,ber_awgn_qpsk2, 'b');


plot(EbNoVec,ber_awgn_bpsk, 'rx:');
xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate');
grid on; drawnow;

title('Performance of BPSK');
legend('Theorical coherent BPSK','Theorical differential BPSK','Theorical coherent QPSK','Theorical differential QPSK','Theorical (BOOK) BPSK');


% BPSK uses two phases (0, 180) to modulate.
% QPSK uses four phases (0, 90, 180, 270) to modulate.
% The BER for BPSK and QPSK are exactly the same
% However, with two bits per symbol in QPSK, the symbol error rate is increased
