function [outData] = nanhandle(data, sampFreq, fmin)

rolloff = 0.95; %Roll-off in percentage
% fmin = 30; %Highpass cutoff

%Find indices of NaNs
nanidxs = find(isnan(data));

if isempty(nanidxs)
    data = data.*tukeywin(length(data), rolloff)';
    outData = highpass(data, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    return;
end

idxs = find(~isnan(data));

%Get all chunk start and end indices
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

%Tukey Window and Highpass filter all chunks
windowedhpchunks = cell(1,length(chunk_end_idxs));
for j = 1:length(chunk_start_idxs)
    tempchunk = data(chunk_start_idxs(j):chunk_end_idxs(j));
    tempchunk = tempchunk.*tukeywin(length(tempchunk), rolloff)';
    tempwindowedchunk = tempchunk;
    tempwindowedhpchunk = highpass(tempwindowedchunk, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    windowedhpchunks{j} = tempwindowedhpchunk;
end

%Estimate the PSD from one of these chunks


end

