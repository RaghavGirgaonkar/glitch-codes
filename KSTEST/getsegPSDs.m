function []= getsegPSDs(launcherparamfilename)


%Load Files
fname = launcherparamfilename;
str = fileread(fname);
launcherparams = jsondecode(str);
hdf5filenames = launcherparams.hdf5filenames;
segfile = launcherparams.segfile;
psdfile = launcherparams.psdfile;
sampFreq = launcherparams.sampFreq;
outfilename = launcherparams.outfilename;
threshold = launcherparams.threshold;
nsegs = 92;

segPSDs = {};
trainidxs = {};
for segnum = 1:nsegs
    %Load Data vector from file(s)
    %Get Segment Data
    segments = textread(segfile, '%s', 'delimiter', '\n');
    Seg = split(segments{segnum});
    segstart = str2num(Seg{2});
    segend = str2num(Seg{3});
    filenums = str2num(Seg{4});
    
    %Get File data and load segment
    if length(filenums) == 1
        files = textread(hdf5filenames, '%s', 'delimiter', '\n');
        fileStr = split(files{filenums});
        filename = fileStr{1};
        GPSstart = str2num(fileStr{2});
        
        %Get segment indices
        startidx = sampFreq*(segstart- GPSstart);
        endidx = sampFreq*(segend- GPSstart);
    
        %Load data from file and create segment data vector
        filedata = h5read(filename, '/strain/Strain')';
        segment = filedata(startidx+1:endidx);
        psdnum = filenums;
    else
        segment = [];
        psdnum = filenums(1);
        for i = 1:length(filenums)
            files = textread(hdf5filenames, '%s', 'delimiter', '\n');
            fileStr = split(files{filenums(i)});
            filename = fileStr{1};
            GPSstart = str2num(fileStr{2});
            GPSend = str2num(fileStr{3});
    
            if GPSend > segstart && GPSstart < segstart && segend > GPSend
                %Get segment indices
                startidx = sampFreq*(segstart- GPSstart);
                endidx = sampFreq*(GPSend -GPSstart);
                %Load data from file and create segment data vector
                filedata = h5read(filename, '/strain/Strain')';
                seg = filedata(startidx:endidx);
                segment = [segment, seg];
            end
            if GPSstart < segend && GPSend > segend && segstart < GPSstart
                %Get segment indices
                startidx = 1;
                endidx = sampFreq*(segend - GPSstart);
                %Load data from file and create segment data vector
                filedata = h5read(filename, '/strain/Strain')';
                seg = filedata(startidx:endidx-1);
                segment = [segment, seg];
            end
        end
    end
    
    Seglen = length(segment)/sampFreq;
    
    %Load PSD and condition data
    %Load PSD
    % psdVecs = load(psdfile);
    % PSD = psdVecs.psd(psdnum,:);
    
    %Condition Data
    rolloff = 0.1;
    fmin = 30;
    % paddedsegment = [zeros(1,sampFreq*rolloff), segment, zeros(1,sampFreq*rolloff)];
    %Window and Highpass filter
    segwin = segment.*tukeywin(length(segment),rolloff*sampFreq/length(segment))';
    seghpass = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    filtsegment = seghpass;
    
    %Find Stationary region for training segment (32s long)
    startidx = 1; 
    seglen = 32*sampFreq;
    endidx = startidx + seglen -1;
    nblocks = 8;
    
    while startidx < length(filtsegment)
        trainsegment = filtsegment(startidx:endidx);
        [hmatrix, ~, ~, ~]=spectogramkstest(trainsegment,sampFreq, nblocks,threshold);
        if hmatrix == 0
            trainstartidx = startidx;
            trainendidx = endidx;
            break;
        else
            startidx = startidx + seglen;
            endidx = startidx + seglen -1;
        end
    end
    
    %Find Welch Estimate
    trainseg = filtsegment(trainstartidx:trainendidx);
    [pxx,f]=pwelch(trainseg, tukeywin(0.5*sampFreq),[],[],sampFreq);
    [interPSD, fvec] = createPSD(sampFreq, Seglen, pxx', f');
    PSD = interPSD./2; %Make two-sided
    
    segPSDs = [segPSDs; PSD];
    trainidxs = [trainidxs; [trainstartidx, trainendidx]];
    
end
save(outfilename,'segPSDs', 'trainidxs');
end
