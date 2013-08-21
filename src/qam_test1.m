%
%   File:      qam_test1.m
%   Author(s): Horacio Sanson
%
EbNoVec = [0:20];

ebno = 10.^(EbNoVec./10);   % Convertir de DB a ratio ordinario

% MATLAB BER for 16QAM
awgn_16qam_matlab = berawgn(EbNoVec, 'qam', 16);
% Equation BER for 16QAM
awgn_16qam_equ = (3/8)*erfc(sqrt((2/5).*ebno)) - (9/64) * (erfc(sqrt((2/5).*ebno)).^2);

% MATLAB BER for 64QAM
awgn_64qam_matlab = berawgn(EbNoVec, 'qam', 64);
% Equation BER for 64QAM
awgn_64qam_equ = (7/24).*erfc(sqrt((1/7).*ebno)) - (49/384).*(erfc(sqrt((1/7).*ebno))).^2;

% Another approximation for 64QAM
%awgn_64qam_equ = (7/24).*erfc(sqrt((1/7).*ebno));

figure(1);
semilogy(EbNoVec,awgn_16qam_matlab, 'm');
hold on;
semilogy(EbNoVec,awgn_16qam_equ, 'm:*');

semilogy(EbNoVec,awgn_64qam_matlab, 'r');
semilogy(EbNoVec,awgn_64qam_equ, 'r:*');
xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate');
grid on; drawnow;

title('Performance of QAM');


