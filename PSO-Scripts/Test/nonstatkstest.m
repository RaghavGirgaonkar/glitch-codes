function [segstart,segend] = nonstatkstest(dataVec,indices)

%Get start and end indices
startidxs = indices(:,1);
endidxs = indices(:,2);

%Get the real part of FFT of each segment
reFFTs = {};
% absFFTs = [];
% threshold = 1.5;
for i  = 1:length(startidxs)
    seg = dataVec(startidxs(i):endidxs(i));
    segFFT = fft(seg);
    realFFT = abs(segFFT);
%     stdFFT = std(realFFT);    
%     rmoutFFT = realFFT((-threshold*stdFFT <= realFFT) & (realFFT <= threshold*stdFFT));
    reFFTs = [reFFTs; realFFT];
%     absFFTs = [absFFTs; abs(segFFT)];
end

%Reject Outliers
realVals = {};
for j= 1:length(startidxs)
    reVals = reFFTs{j};
    absreVals = abs(reVals);
    [N,edges] = histcounts(absreVals);
    CSM = cumsum(N);
    idx = find(CSM<=0.995*length(absreVals), 1, 'last');
    threshold = edges(idx);
    clndreVals = reVals((-threshold <= reVals) & (reVals <= threshold));
    realVals = [realVals; clndreVals];
    disp(size(clndreVals));
%     disp(num2str(j));
end
% realVals = reFFTs;

%Form KS test matrix
pmatrix = zeros(length(endidxs));
hmatrix = zeros(length(endidxs));
kmatrix = zeros(length(endidxs));
for i = 1:length(startidxs)
    for j = i:length(startidxs)
        [h,p,k] = kstest2(realVals{i},realVals{j});
        pmatrix(i,j) = p;
        pmatrix(j,i) = p;
        hmatrix(i,j) = h;
        hmatrix(j,i) = h;
        kmatrix(i,j) = k;
        kmatrix(j,i) = k;
    end
end

end

