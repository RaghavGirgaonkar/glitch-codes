function [whtndfiltdata] = correctwhtndstrain(whtndfiltdata, sampFreq,nanchunk_start_idxs, nanchunk_end_idxs)
    N = length(whtndfiltdata);
    n = 0.6;
    for k = 1:length(nanchunk_start_idxs)
        whtndfiltdata(nanchunk_start_idxs(k):nanchunk_end_idxs(k)) = randn(1,nanchunk_end_idxs(k) - nanchunk_start_idxs(k)+1);
    end
    for i = 1:length(nanchunk_end_idxs)
        if nanchunk_start_idxs(i) ~= 1 && nanchunk_start_idxs(i)/sampFreq >= n
            nanidx = nanchunk_start_idxs(i);
            whtndfiltdata(nanidx - floor(n*sampFreq): nanidx + floor(n*sampFreq)) = randn(1, floor(2*n*sampFreq));
        end
        if nanchunk_end_idxs(i) ~= N && (N - nanchunk_end_idxs(i))/sampFreq >= n
            nanidx = nanchunk_end_idxs(i);
            whtndfiltdata(nanidx - floor(n*sampFreq): nanidx + floor(n*sampFreq)) = randn(1, floor(2*n*sampFreq));
        end
    end
end

