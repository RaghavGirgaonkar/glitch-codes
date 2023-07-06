function [segstart,segend] = nonstatkstest(dataVec,indices)

%Get start and end indices
startidxs = indices(:,1);
endidxs = indices(:,2);

%Get the real part of FFT of each segment
reFFTs = {};
% absFFTs = [];
threshold = 1.5;
for i  = 1:length(startidxs)
    seg = dataVec(startidxs(i):endidxs(i));
    segFFT = fft(seg);
    realFFT = real(segFFT);
    stdFFT = std(realFFT);    
    rmoutFFT = realFFT((-threshold*stdFFT <= realFFT) & (realFFT <= threshold*stdFFT));
    reFFTs = [reFFTs; rmoutFFT];
%     absFFTs = [absFFTs; abs(segFFT)];
end

%Reject Outliers outside of ~ 1.5 sigma
% stdFFT = std(reFFTs);
% threshold = 1.5;
% rmoutFFT = reFFTs((-threshold*stdFFT <= reFFTs) & (reFFTs <= 1.5*threshold));


%Form KS test matrix
pmatrix = zeros(length(endidxs));
hmatrix = zeros(length(endidxs));
kmatrix = zeros(length(endidxs));
for i = 1:length(startidxs)
    for j = i:length(startidxs)
        [h,p,k] = kstest2(reFFTs{i},reFFTs{j});
        pmatrix(i,j) = p;
        pmatrix(j,i) = p;
        hmatrix(i,j) = h;
        hmatrix(j,i) = h;
        kmatrix(i,j) = k;
        kmatrix(j,i) = k;
    end
end

end

