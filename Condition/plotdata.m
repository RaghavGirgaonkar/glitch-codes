function []=plotdata(outData, nanchunk_start_idxs, nanchunk_end_idxs, chunk_start_idxs, chunk_end_idxs)
    nanseries = zeros(size(outData));
    dataseries = zeros(size(outData));

    for k = 1:length(nanchunk_start_idxs)
        nanseries(nanchunk_start_idxs(k):nanchunk_end_idxs(k)) = outData(nanchunk_start_idxs(k):nanchunk_end_idxs(k));
    end

    for j = 1:length(chunk_start_idxs)
        dataseries(chunk_start_idxs(j):chunk_end_idxs(j)) = outData(chunk_start_idxs(j):chunk_end_idxs(j));
    end

    figure;
    plot(dataseries, Color='blue', DisplayName='Original');
    hold on;
    plot(nanseries, Color='red', DisplayName='NaN-filled regions');
    legend;
    hold off;
end