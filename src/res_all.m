%
%   File:      res_all.m
%   Author(s): Horacio Sanson
%   Revision : 2007/11/13

%   Description:
%     Wrapper script that executes all other res_scrXX.m scripts.

%% Probability of resolution vs SNR for different filter sizes
%  res_scd01       % WOSA
%  res_scd02       % MVDR
%  res_scd03       % NMVDR
%  res_scd04       % APES


%% Probability of resolution vs SNR for different methods
%  res_scd06       % df = 1
%  res_scd08       % df = 0.5
%  res_scd10       % df = 0.25


%% Probability of resolution vs frequency separation
res_scd21       % WOSA
res_scd22       % MVDR
res_scd23       % NMVDR
res_scd24       % APES

%% Probability of resolution for different frequency separations 
% SNR = 20 dB
res_scd26
