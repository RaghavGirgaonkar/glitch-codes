function []=datacond(filename, outfilename,tstart, seglen, sampFreq, winlen)

fmin = 30;
% sampfreq = 4096;

data_file = h5read(filename, '/strain/Strain')';

%Account for NaNs in data, find all non-NaN indices
nanidxs = ~isnan(data_file);

idxs = find(~isnan(data_file));

data = data_file(nanidxs);

trainingidxs = tstart*sampFreq:(tstart+seglen)*sampFreq;

if ~ismember(trainingidxs,idxs)
   error("Choose different segment for Welch estimate");
end


N = length(data);
Tsig = N/sampFreq;
timeVec = (0:N-1)*(1/sampFreq);

%Highpass Filter the data
% [b,a] = butter(8,fmin/(sampFreq/2),'high');

%Filter Data

filtdata = highpass(data, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);

% filtdata = filtdata_temp(t*sampFreq+1:end - t*sampFreq);

%Take welch estimate of specified segment
if ~isempty(idx)
    tempdata = filtdata((tstart - nantime)*sampFreq: (tstart- nantime+seglen)*sampFreq);
else
    tempdata = filtdata(tstart*sampFreq: (tstart+seglen)*sampFreq);
end

%Zero-pad segment before Welch estimate
% tempdata_t = [zeros(1,winlen*sampFreq), tempdata, zeros(1,winlen*sampFreq)];

[pxx,f]=pwelch(tempdata, tukeywin(winlen*sampFreq),[],[],sampFreq);

logwelchPSD = log10(pxx');

%Interpolate Welch Estimate to create entire PSD vector
PSD = createPSD(sampFreq, Tsig, logwelchPSD, f');

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

createHDF5file(outfilename, whtndfiltdata, TF, N);