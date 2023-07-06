function [] = getpsds(hdf5filenames, psdfilename, sampFreq, trainlen, seglen)

%Get Filenames
files = textread(hdf5filenames, '%s', 'delimiter', '\n');
psd = [];
fmin = 30;
for i = 1:length(files)
    fileStr = split(files{i});
    filename = fileStr{1};
    filedata = h5read(filename, '/strain/Strain')';
    filelen = length(filedata)/sampFreq;

    %Get first 64 seconds to estimate PSD
    traindata = filedata(1:64*sampFreq);

    %Highpass 
    traindatahpass = highpass(traindata, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);

    %Pwelch
    [pxx, f] = pwelch(traindatahpass, tukeywin(trainlen*sampFreq), [], [], sampFreq);
    
%     Interpolate to match file lengths
    logwelchPSD = log10(pxx');
    [PSD, ~] = createPSD(sampFreq, seglen, logwelchPSD, f');
    
    %Add to Matrix
    psd = [psd;PSD];

end

%Save PSD file
save(psdfilename, 'psd');

end

