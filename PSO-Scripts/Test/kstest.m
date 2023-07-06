function [segstart,segend] = kstest(dataVec,indices)

%Get start and end indices
startidxs = indices(:,1);
endidxs = indices(:,2);

%Get the real part of FFT of each segment
reFFTs = [];
absFFTs = [];
for i  = 1:length(startidxs)
    seg = dataVec(startidxs(i):endidxs(i));
    segFFT = fft(seg);
    reFFTs = [reFFTs; real(segFFT)];
    absFFTs = [absFFTs; abs(segFFT)];
end

%Reject Outliers outside of ~ 1.5 sigma
stdFFT = std(reFFTs);
threshold = 10;
rmoutFFT = reFFTs((-threshold*stdFFT <= reFFTs) & (reFFTs <= 1.5*threshold));


%Form KS test matrix
pmatrix = zeros(length(endidxs));
hmatrix = zeros(length(endidxs));
kmatrix = zeros(length(endidxs));
for i = 1:length(startidxs)
    for j = i:length(startidxs)
        [h,p,k] = kstest2(rmoutFFT(i,:),rmoutFFT(j,:));
        pmatrix(i,j) = p;
        pmatrix(j,i) = p;
        hmatrix(i,j) = h;
        hmatrix(j,i) = h;
        kmatrix(i,j) = k;
        kmatrix(j,i) = k;
    end
end

end

