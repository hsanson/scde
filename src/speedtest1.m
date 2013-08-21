%
%   File:      speedtest1.m
%   Author(s): Horacio Sanson

%   Description:
%     Speed test of several matrix inversion methods
%     Speed test of Capon spectral estimation methods

clear;clc;

x = randn(1024,1);      % random signal
R = corr_matlab(x,x,length(x)/2);

%% Inversion using matlab inv method
tic;
Ri = inv(R);
inv_time = toc

%% Inversion using matlab \ operator
tic;
Ri = R\eye(size(R));
back_time = toc

%% Inverion using lu

tic;
[L U] = lu(R);
Ri2 = (U\eye(size(U)))*(L\eye(size(L)));
LU_time = toc

%% Inverion using cholesky
tic;
L = chol(R);
Ri3 = (L\eye(size(L)))*(L'\eye(size(L)));
cholesky_time = toc


%% Capon methods speed measure

tic;
C1 = capon(x,x,length(x)/8);
capon_time = toc

tic;
C2 = capon2(x,x,length(x)/8);
capon2_time = toc

tic;
C3 = caponL(x,x,length(x)/8);
caponL_time = toc

tic;
C3 = caponL2(x,x,length(x)/8);
caponL2_time = toc

tic;
C3 = apes(x,x,20);
apes_time = toc


%% scd speed test 
tic;
S = scd(x,512,'WOSA');
scd_wosa = toc

tic;
S = scd(x,512,'MVDR');
scd_mvdr = toc

tic;
S = scd(x,512,'NMVDR');
scd_nmvdr = toc

tic;
S = scd(x,512,'APES');
scd_apes = toc

