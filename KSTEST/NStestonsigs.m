function NStestonsigs(filtdata, psd, trainidxs, sampFreq, snr, threshold, outfilename)

% %Read file/data vector
% data = h5read(filename, '/strain/Strain')';
% 
% %Highpass data vector about fmin = 30 Hz
% rolloff = 4;
% fmin = 30; sampFreq = 4096;
% segment = data;
% 
% %Window and Highpass filter
% segwin = segment.*tukeywin(length(segment),rolloff*sampFreq/length(segment))';
% seghpass = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
% filtdata = seghpass;
% 
% %Estimate PSD and create final PSD vector
% trainidxs = [1,64*sampFreq];
% [pxx_1,f_1] = pwelch(filtdata(1:64*sampFreq), tukeywin(4*sampFreq),[],[],sampFreq);
% [psd,~] = createPSD(sampFreq, 4096, pxx_1', f_1');

%Create signal to inject with user specified SNR
frange = [30,700]; 
masses = [2,2]; 
ta = 138; 
siglen = 512; 
initial_phase = 0; 
phase = 0;
T_sig = 54;

[signal] = snrsiginj(psd, snr, masses, frange, ta, phase, siglen, sampFreq, initial_phase);

[~, whtndstd, ~] = segdatacond(filtdata, psd, sampFreq, trainidxs);

sigstrain = whtndstd*signal;

%Add to data strain
filtdata = filtdata + sigstrain;

%Use injected signal data segment to estimate PSD 
trainseg = filtdata((ta + 15)*sampFreq:(ta + 15 + 32)*sampFreq);

[pxx,f] = pwelch(trainseg, tukeywin(2*sampFreq),[],[],sampFreq);
[PSD,~] = createPSD(sampFreq, siglen, pxx', f');

%% Run spectogram based KSTest to see if injected signal is detected
[hmatrix, filtlogS, pmatrix, kmatrix]=spectogramkstest(trainseg,sampFreq, 8,threshold);

if hmatrix == 0
    detected = 0;
else
    detected = 1;
end

% Get maximum KS distance value and minimum pvalue

maxkdist = max(kmatrix{1}(:));

minpvalue = min(pmatrix{1}(:));

%Create final whitened strain

[whtndseg,~, TFtotal] = segdatacond(filtdata, PSD, sampFreq, trainidxs);

%% Retrieve signal from Matched-Filtering
fmin = frange(1);
fmax = frange(2);
N = siglen*sampFreq;
fpos = (0:floor(N/2))*(1/siglen);
[A,avec, phaseDiff] = preprocessing(fmin,fmax, fpos, siglen, N);

dataY = whtndseg;

AbysqrtPSD = A.*TFtotal;

%% Create General Normalization Factor
% Scalar factor of 1/N is due to Parseval's theorem
dataLen = N;
innProd = (1/dataLen)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
genNormfac = 1/sqrt(genNormfacSqr);

%% Data Products 
fftdataY = fft(dataY);
fftdataY = fftdataY.*A;%Pre-multiply frequency magnitude vector A for optimization

%% Get FFT of data by total PSD
fftdataYbyPSD = fftdataY.*TFtotal;

%% Run matched filtering
m1 = masses(1);
m2 = masses(2);
phaseq0 = gen2PNwaveform(fpos, 0, 0, 30, 700, m1,...
2,siglen,0,1,N,avec, genNormfac);
fftq0 = phaseq0;
fftq1 = phaseq0.*phaseDiff;
mf1 = matchedfiltering(fftdataYbyPSD, fftq0);
mf2 = matchedfiltering(fftdataYbyPSD, fftq1);

mftimeseries = sqrt(mf1(1:end - T_sig*sampFreq).^2 + mf2(1:end - T_sig*sampFreq).^2);

%Get SNR

[maxval, maxarg] = max(mftimeseries);

% mfstd = std(mftimeseries(1:100000));

% retSNR = maxval/(mfstd*sqrt(2));

retSNR = maxval;


%Write results to file
results = [snr, retSNR, detected, maxkdist, minpvalue];

fileID = fopen(outfilename,'a');
fprintf(fileID,'%f\t%f\t%d\t%f\t%f\n',results);
fclose(fileID);

end


















