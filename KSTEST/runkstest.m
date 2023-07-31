function [Hfinal,filespctgrm] = runkstest(filename, sampFreq)
%RUNKSTEST 

%Load file 
data = h5read(filename, '/strain/Strain')';

%Highpass data 
rolloff = 1;
fmin = 30; 
%Window and Highpass filter
segwin = data.*tukeywin(length(data),rolloff*sampFreq/length(data))';
filtdata = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);

%Run 2 sample KS-test on non-overlapping blocks of 32 second duration
Hfinal = [];
filespctgrm = [];
% datalen = length(data)/sampFreq;
startidx = 1; 
seglen = 32*sampFreq;
endidx = startidx + seglen -1;
nblocks = 8;
threshold = 0.13;

while startidx < length(data)
    segment = filtdata(startidx:endidx);
%     segment = data(startidx:endidx);
%     %Window and Highpass segment
%     segwin = segment.*tukeywin(length(segment),rolloff*sampFreq/length(segment))';
%     filtsegment = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    [hfinal, outspctgrm, ~, ~] = spectogramkstest(segment, sampFreq, nblocks, threshold);

    %Get Hmatrix time info
    htime = zeros(1,length(hfinal));
    for i = 1:length(hfinal)
        htemp = htime;
        htemp(i:end) = hfinal(i,i:end);
        htime = htemp | htime;
    end
    
    %Convert from 4 second blocks to 1 second
    HTIME = [];
    for i = 1:length(htime)
        for j = 1:4
            HTIME = [HTIME, htime(i)];
        end
    end

    Hfinal = [Hfinal, HTIME];
    filespctgrm = [filespctgrm, outspctgrm];
    startidx = startidx + seglen;
    endidx = startidx + seglen -1;
end


end

