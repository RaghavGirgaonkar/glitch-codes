function [PSD, fvec]=createPSD(sampFreq, Tsig, logwelchPSD, freqs)
%% Script to create interpolated PSD from SHAPES and PWELCH Estimates

%% Data Parameters
Fs = sampFreq;
T = Tsig;
N = Fs*T;

kNyq = floor(N/2);
fvec = (0:(kNyq))*Fs/N;

% nyqfreq = winlen*Fs;

% freqs = (0:nyqfreq)*(Fs/(2*nyqfreq));

%% 1-D Interpolation

loginterPSD = interp1(freqs, logwelchPSD, fvec);

% %% Antilog

PSD = (10.^loginterPSD)./2;
% PSD = (loginterPSD)./2;
