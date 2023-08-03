function [hfinal, filtlogS, Pmatrixs, Kmatrixs]=spectogramkstest(dataVec,sampFreq, nblocks, threshold)

    %Make spectrogram first
    [s,f,~] = spectrogram(dataVec, 1024, [],[], sampFreq);
    logS = log10(abs(s));
    [~, lenT] = size(logS);

    %Get blocks corresponding to freqs between 30 and 700 Hz
    freqs = find(50<=f & f <=700)';
    filtlogS = logS(freqs,:);
    
    %Remove lines
    filtlogS = [filtlogS(1:62,:);filtlogS(69:110,:);filtlogS(119:end,:)];
    
    %Normalize filtered spectogram
    avgVal = median(filtlogS(:));
%     stdVal = std(filtlogS(:));
%     filtlogS = (filtlogS - avgVal)./stdVal;
    filtlogS(filtlogS < avgVal) = avgVal;
    maxVal = max(filtlogS(:));
    filtlogS = filtlogS/(-1*maxVal);

     

    jump = floor(lenT/nblocks);

    %Make Fragments and take mean along time axis for each frag
    Hmatrixs = {};
    Pmatrixs = {};
    Kmatrixs = {};
    for l = 1:mod(size(filtlogS,2),nblocks)
        specavgs = {};
         specstart = l; 
         specend = specstart + jump - 1;
         while specstart < lenT
    %         disp(num2str(specstart));
    %         disp(num2str(specend));
            frag = filtlogS(:, specstart:specend);
            fragavg = mean(frag,2);
    %         fragavg = median(frag,2);
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
    %     hmatrix = zeros(nblocks);
        kmatrix = zeros(nblocks);
        for i = 1:nblocks
            for j = i:nblocks
                [~,p,k] = kstest2(specavgs{i},specavgs{j});
                pmatrix(i,j) = p;
                pmatrix(j,i) = p;
    %             hmatrix(i,j) = h;
    %             hmatrix(j,i) = h;
                kmatrix(i,j) = k;
                kmatrix(j,i) = k;
            end
        end
        hmatrix = kmatrix;
        hmatrix(hmatrix > threshold) = 1; hmatrix(hmatrix < threshold) = 0;
        Hmatrixs = [Hmatrixs; hmatrix];
        Pmatrixs = [Pmatrixs; pmatrix];
        Kmatrixs = [Kmatrixs; kmatrix];
    end
    
    y = 1;
    for i = 1:length(Hmatrixs)
        if Hmatrixs{i} == 0
            y = 0;
        end
    end
    if y == 0
        hfinal = zeros(nblocks);
    else
        hfinal = zeros(nblocks);
        for i = 1: length(Hmatrixs)
            hfinal = hfinal | Hmatrixs{i};
        end
    end
end

