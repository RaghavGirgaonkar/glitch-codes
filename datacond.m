function []=datacond(filename, outfilename,tstart, seglen, sampFreq, winlen)

fmin = 30;
% sampfreq = 4096;

data = h5read(filename, '/strain/Strain')';

N = length(data);
Tsig = N/sampFreq;
timeVec = (0:N-1)*(1/sampFreq);

%Highpass Filter the data and handle NaN-regions
[filtdata, PSD, nanchunk_start_idxs, nanchunk_end_idxs] = nanhandle(data, sampFreq, fmin, tstart, seglen, winlen);

%Create Transfer Function
TF = 1./sqrt(PSD);

%Set Values of TF below fmin to 0
TF(1:fmin*sampFreq) = 0;

%Create entire Transfer Function vector
negFStrt = 1-mod(N,2);
kNyq = floor(N/2)+1;

TFtotal = [TF, TF((kNyq-negFStrt):-1:2)];

%Whiten Filter Data with Transfer Function
%Window data before FFT with a Tukey-window with a 0.5 seconds rolloff

rolloff = 4; %Roll-off in seconds
winfiltdata = filtdata.*tukeywin(length(filtdata), rolloff*sampFreq/N)';

fftfiltdata = fft(winfiltdata);

whtndfftfiltdata = fftfiltdata.*TFtotal;

whtndfiltdata = ifft(whtndfftfiltdata);

%Divide by variance of whitened strain so that final whitened vector has
%unit strain
whtndstd = std(whtndfiltdata(tstart*sampFreq: (tstart+seglen)*sampFreq));
whtndfiltdata = whtndfiltdata/whtndstd;

%Correct edge effects for files with NaNs
if ~isempty(isnan(data))
    whtndfiltdata = correctwhtndstrain(whtndfiltdata, sampFreq, nanchunk_start_idxs, nanchunk_end_idxs);
end
createHDF5file(outfilename, whtndfiltdata, TF, N);