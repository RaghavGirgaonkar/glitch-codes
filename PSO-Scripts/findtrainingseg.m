function [segStart, segEnd]= findtrainingseg(filename, fmin, sampFreq, seglen, overlap)
    %Load HDF5 file as row vector
    data = h5read(filename,'/strain/Strain')';
    %Step 1: Tukey-Window and Highpass data above <fmin> Hz
    rolloff = 1;
    datalen = length(data)/sampFreq;
    datawin = data.*tukeywin(length(data), rolloff*sampFreq/datalen)';
    datawinhpass = highpass(datawin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
    
    %Step 2: Divide Highpassed data into overlapping segments and run KS
    %test on the FFT of each segment
    segStart = 1;
    segEnd = seglen*sampFreq;
    segSamples = seglen*sampFreq;
    %Get seg K FFT
    segrolloff = 0.5;
    segK = datawinhpass(segStart:segEnd);
    segKwin = segK.*tukeywin(length(segK),segrolloff*sampFreq/seglen)';
    segKFFT = abs(fft(segKwin));
    segKFFT = segKFFT(1:floor(segSamples/2));
    %Get next 5 overlapping segments and their FFT
    n = 5 ;
    overlapsegFFT = zeros(n,floor(segSamples/2));
    for i = 1:n
        Sstart = segEnd - overlap*sampFreq;
        Send = Sstart + seglen*sampFreq-1;
        SegI = datawinhpass(Sstart:Send);
        SegIwin = SegI.*tukeywin(length(SegI),segrolloff*sampFreq/seglen)';
        SegIFFT = abs(fft(SegIwin));
        SegIFFT = SegIFFT(1:floor(segSamples/2));
        overlapsegFFT(i,:) = SegIFFT;
    end
    %KS Test with SegK and SegIs
    kvals = zeros(1,n);
    hvals = zeros(1,n);
    pvals = zeros(1,n);
    for j = 1:n
        [h,p,k] = kstest2(segKFFT, overlapsegFFT(j,:));
        kvals(j) = k;
        hvals(j) = h;
        pvals(j) = p;
    end
end