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
segments = textread(segfile, '%s', 'delimiter', '\n');
nsegs = length(segments);

segPSDs = {};
trainidxs = {};
for segnum = 255:258
    disp(segnum);
    %Load Data vector from file(s)
    %Get Segment Data
    segments = textread(segfile, '%s', 'delimiter', '\n');
    Seg = split(segments{segnum});
    segstart = str2num(Seg{2});
    segend = str2num(Seg{3});
    filenums = str2num(Seg{4});
%     disp(segstart)
%     disp(segend)
    
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
    O1timeline = [1126051217, 1137254417];
    O2timeline = [1164556817, 1187733618];
    O3atimeline = [1238166018, 1253977218];
    O3btimeline = [1256655618, 1269363618];

    startidx = 1; 
    seglen = 32*sampFreq;
    endidx = startidx + seglen -1;
    nblocks = 8;
    ntrialsthresh = 5;
    O1linefreqs = [30,46; 54,66;118,126;178,182;305,340;498,510;626,640];
    O2linefreqs = [50,70; 298,322;494,502;505,520];
    O3linefreqs = [302,318;500,520];

    

    if segstart >= O1timeline(1) && segstart < O1timeline(2)
        linefreqs = O1linefreqs;
        disp('O1File');
    elseif segstart >= O2timeline(1) && segstart < O2timeline(2)
        linefreqs = O2linefreqs;
        disp('O2File');
    else
        linefreqs = O3linefreqs;
        disp('O3File');
    end

    inParams = struct('sampFreq',sampFreq,...
                       'frange',[30,700],...
                       'linefreqs',linefreqs,...
                       'specres',1024,...
                       'nblocks',nblocks,...
                       'threshold', threshold,...
                       'ntrialsthresh', ntrialsthresh);

    while startidx < length(filtsegment)
        trainsegment = filtsegment(startidx:endidx);
        [hmatrix]=NSPECKT(trainsegment,inParams);
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
