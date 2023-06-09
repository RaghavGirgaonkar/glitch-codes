function outNoise = genColNoise(nSamples,psdVals,fltrOrdr,sampFreq)
    %Generate a realization of stationary Gaussian noise with given 2-sided PSD
%Y = genColNoise(N,PSD,Ordr,Fs)
%Generates a realization Y of stationary gaussian noise with a target
%2-sided power spectral density given by PSD. Fs is the sampling frequency
%of Y. PSD is an M-by-2 matrix containing frequencies and the corresponding
%PSD values in the first and second columns respectively. The frequencies
%must start from 0 and end at Fs/2. The order of the FIR filter to be used
%is given by Ordr.

%Soumya D. Mohanty, Mar 2019

% Added Tukey-windowing option 

% Raghav Girgaonkar, Apr 2023

%% Design FIR filter with T(f)= square root of target PSD
freqVec = psdVals(:,1);
sqrtPSD = sqrt(psdVals(:,2));
b = fir2(fltrOrdr,freqVec/(sampFreq/2),sqrtPSD,tukeywin(fltrOrdr+1));


%% Uncomment following block of code if creating noise vector
%% leave commented otherwise
% (Comment out the line below if new realizations of WGN are needed in each run of this script)
% rng('default'); 

inNoise = randn(1,nSamples);
% inNoise = [randn(1,fltrOrdr), inNoise_t, randn(1,fltrOrdr)];
outNoise = sqrt(sampFreq)*fftfilt(b,inNoise);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

