function [] = makesegments(hdf5filenames, seglen, overlap)

%Get number of segments per file
totalLen = 4096;

nSeg = floor((totalLen-seglen)/(seglen - overlap));

%Get filenames
files = textread(hdf5filenames, '%s', 'delimiter', '\n');
outputFile = fullfile('segments.txt');
fid = fopen(outputFile, 'w');
counter = 1;
fileStr = split(files{1});
GPSstart = str2num(fileStr{2});
GPSend = str2num(fileStr{3});
% name = fileStr{1};
segStart = GPSstart;
segEnd = segStart + seglen;
for i = 1:length(files)
    fileStr = split(files{i});
    GPSstart = str2num(fileStr{2});
    GPSend = str2num(fileStr{3});

    %Make segments
    filecounts = [i];
    while segEnd < GPSend
        fprintf(fid, '%d\t%d\t%d\t[%s]\n', counter, segStart, segEnd, join(string(filecounts),','));
        segStart = segEnd - overlap;
        segEnd = segStart + seglen;
        counter = counter + 1;
    end
    
    if segEnd ~= GPSend && i == length(files)
        fprintf(fid, '%d\t%d\t%d\t[%s]\n', counter, segStart, GPSend, join(string(filecounts),','));
        return;
    end

    while i < length(files) && segStart < GPSend && segEnd >= GPSend
        if segEnd == GPSend
            fprintf(fid, '%d\t%d\t%d\t[%s]\n', counter, segStart, segEnd, join(string(filecounts),','));
            segStart = segEnd - overlap;
            segEnd = segStart + seglen;
            counter = counter + 1;
            continue;
        else
            filecounts = [i,i+1];
            fprintf(fid, '%d\t%d\t%d\t[%s]\n', counter, segStart, segEnd, join(string(filecounts),','));
            segStart = segEnd - overlap;
            segEnd = segStart + seglen;
            counter = counter + 1;
            continue;
        end

    end
    
    if segStart >= GPSend && segEnd >= GPSend
        continue;
    end

end
fclose(fid);

end

