function [hfinal, Hmatrixs, Pmatrixs, Kmatrixs]=spectogramkstest2(dataVec,sampFreq, nblocks)

    %Make spectrogram first
    [s,f,~] = spectrogram(dataVec, 1024, [],[], sampFreq);
    logS = log10(abs(s));
    [~, lenT] = size(logS);

    %Normalize Spectogram between 30 and 700 Hz
    bandfreqs = find(30<=f & f <=700)';
    bandlogS = logS(bandfreqs,:);
    avgVal = median(bandlogS(:));
    logS(logS < avgVal) = avgVal;
    maxVal = max(bandlogS(:));
    logS = logS/(-1*maxVal);

    %Iterate over small bandwidths in the frequency range
    Hmatrixs = {};
    Pmatrixs = {};
    Kmatrixs = {};
    bandwidth = 250; %Bandwidth of subblock in Hz
    freqjump = 250; %Jump from one block to next in frequency domain in Hz
    fmin = 50; fmax = fmin + bandwidth;
    
    fFinal = 700;

    while fmin < fFinal
        
%         freqrange = [fmin, fmax];
%         disp(freqrange);
        %Get blocks corresponding to freqs between 30 and 700 Hz
        freqs = find(fmin<=f & f <=fmax)';
        filtlogS = logS(freqs,:);
        
        %Normalize filtered spectogram
%         avgVal = mean(filtlogS(:));
%         filtlogS(filtlogS < avgVal) = avgVal;
%         maxVal = max(filtlogS(:));
%         filtlogS = filtlogS/(-1*maxVal);
    
        jump = floor(lenT/nblocks);
    
        %Make Fragments and take mean along time axis for each frag
        specavgs = {};
         specstart = 1; 
         specend = specstart + jump - 1;
         while specstart < lenT 
    %         disp(num2str(specstart));
    %         disp(num2str(specend));
            frag = filtlogS(:, specstart:specend);
            fragavg = mean(frag,2);
%             fragavg = median(frag,2);
            specavgs = [specavgs;fragavg];
            specstart = specend;
            if specstart + jump > lenT
                specend = lenT;
            else
                specend = specstart + jump - 1;
            end
         end
    
         %Run KSTest
        pmatrix = zeros(nblocks);
        hmatrix = zeros(nblocks);
        kmatrix = zeros(nblocks);
        for i = 1:nblocks
            for j = i:nblocks
                [h,p,k] = kstest2(specavgs{i},specavgs{j}, 'Alpha', 0.01);
                pmatrix(i,j) = p;
                pmatrix(j,i) = p;
                hmatrix(i,j) = h;
                hmatrix(j,i) = h;
                kmatrix(i,j) = k;
                kmatrix(j,i) = k;
            end
        end
        %Find and display scatter of p values
        pvals = [];
        for i = 1:nblocks
            for j = i+1:nblocks
                pvals = [pvals, pmatrix(i,j)];
            end
        end
        disp([fmin, fmax,std(pvals)])
        %Append h p and k matrices to corresponding cell arrays
        Hmatrixs = [Hmatrixs; hmatrix]; 
        Pmatrixs = [Pmatrixs; pmatrix];
        Kmatrixs = [Kmatrixs; kmatrix];
        fmin = fmin + freqjump;
        fmax = fmax + freqjump;
        if fmin + freqjump > fFinal
            fmax = fFinal;
        end
    end

    %Calculate final Hmatrix by OR operations on all previous h matrices
    hfinal = zeros(nblocks);
    for i = 1: length(Hmatrixs)
        hfinal = hfinal | Hmatrixs{i};
    end
end

