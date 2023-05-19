function [outData, PSD] = nanhandle(data, sampFreq, fmin, tstart, seglen, winlen)

rolloff = 0.99; %Roll-off in percentage
% fmin = 30; %Highpass cutoff
N = length(data);
Tsig = N/sampFreq;
timeVec = (0:N-1)*(1/sampFreq);

%Find indices of NaNs
nanidxs = find(isnan(data));

if isempty(nanidxs)
%     data = data.*tukeywin(length(data), rolloff)';
    outData = highpass(data, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    tempdata = outData(tstart*sampFreq: (tstart+seglen)*sampFreq);
    [pxx,f]=pwelch(tempdata, tukeywin(winlen*sampFreq),[],[],sampFreq);

    logwelchPSD = log10(pxx');

    %Interpolate Welch Estimate to create entire PSD vector
    [PSD, ~] = createPSD(sampFreq, Tsig, logwelchPSD, f');
    return;
end

idxs = find(~isnan(data));

%% Get all Non-NaN chunk start and end indices
chunk_start_idxs = [];
chunk_end_idxs = [];

for i = 1:length(idxs)-1
    if i == 1
        chunk_start_idxs = [chunk_start_idxs, idxs(i)];
        continue;
    end

    if idxs(i+1) - idxs(i) ~= 1
        chunk_end_idxs = [chunk_end_idxs, idxs(i)];
        continue;
    end

    if idxs(i) - idxs(i-1) ~= 1
        chunk_start_idxs = [chunk_start_idxs, idxs(i)];
        continue;
    end
end

if length(chunk_start_idxs) ~= length(chunk_end_idxs)
   chunk_end_idxs = [chunk_end_idxs, idxs(end)]; 
end

%% Get all NaN chunk start and end indices
nanchunk_start_idxs = [];
nanchunk_end_idxs = [];

for i = 1:length(nanidxs)-1
    if i == 1
        nanchunk_start_idxs = [nanchunk_start_idxs, nanidxs(i)];
        continue;
    end

    if nanidxs(i+1) - nanidxs(i) ~= 1
        nanchunk_end_idxs = [nanchunk_end_idxs, nanidxs(i)];
        continue;
    end

    if nanidxs(i) - nanidxs(i-1) ~= 1
        nanchunk_start_idxs = [nanchunk_start_idxs, nanidxs(i)];
        continue;
    end
end

if length(nanchunk_start_idxs) ~= length(nanchunk_end_idxs)
   nanchunk_end_idxs = [nanchunk_end_idxs, nanidxs(end)]; 
end

%% Tukey Window and Highpass filter all chunks
windowedhpchunks = cell(1,length(chunk_end_idxs));
for j = 1:length(chunk_start_idxs)
    tempchunk = data(chunk_start_idxs(j):chunk_end_idxs(j));
%     tempchunk = tempchunk.*tukeywin(length(tempchunk), rolloff)';
    tempwindowedchunk = tempchunk;
    tempwindowedhpchunk = highpass(tempwindowedchunk, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    windowedhpchunks{j} = tempwindowedhpchunk;
end

%% Estimate the PSD from one of these chunks from user specified segment
trainstartidx = tstart*sampFreq;
trainendidx = (tstart+seglen)*sampFreq;

%Check if training segment is feasible
if sum(ismember(idxs, trainstartidx:trainendidx)) ~= length(trainstartidx:trainendidx)
       error('Training Segment specified contains NaNs, choose different segment');
end

%Find which chunk contains the training segment
chunkidx = find(chunk_start_idxs <= trainstartidx, 1,'last');

%Get training segment and Estimate PSD by PWelch
trainingdata = windowedhpchunks{chunkidx}(trainstartidx - chunk_start_idxs(chunkidx)+1: trainendidx - chunk_start_idxs(chunkidx));

[pxx,f]=pwelch(trainingdata, tukeywin(winlen*sampFreq),[],[],sampFreq);
logwelchPSD = log10(pxx');

%% Create Colored noise to fill in NaN regions

%Get size of NaN regions in seconds
noiselen = length(nanidxs)/sampFreq;

%Interpolate Pwelch to match this length
[PSD_temp, fvec_temp] = createPSD(sampFreq, noiselen, logwelchPSD, f');
[PSD, ~] = createPSD(sampFreq, Tsig, logwelchPSD, f');

%Create Colored noise using this PSD
fltrOrdr = 10000;
n = 4; %Extra data in seconds to add to remove later
outNoise_temp = genColNoise((noiselen+2*n)*sampFreq, [fvec_temp(:), PSD_temp(:)], fltrOrdr, sampFreq);
outNoise = outNoise_temp(n*sampFreq+1: end - n*sampFreq);

%% Create Final Vector with Tukey-Windowed Chunks and NaN regions filled
outData = zeros(size(data));

%First Add Tukey Windowed Chunks
for k = 1:length(chunk_start_idxs)
    outData(chunk_start_idxs(k):chunk_end_idxs(k)) = windowedhpchunks{k};
end

%Lastly add Generated Colored noise in place of NaNs
a = 1;
for k = 1:length(nanchunk_start_idxs)
    outData(nanchunk_start_idxs(k):nanchunk_end_idxs(k)) = outNoise(a:a+nanchunk_end_idxs(k) - nanchunk_start_idxs(k));
    a = a+nanchunk_end_idxs(k) - nanchunk_start_idxs(k) + 1;
end
plotdata(outData, nanchunk_start_idxs, nanchunk_end_idxs, chunk_start_idxs, chunk_end_idxs);
end

