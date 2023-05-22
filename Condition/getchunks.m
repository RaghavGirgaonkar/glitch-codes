function [chunk_start_idxs, chunk_end_idxs] = getchunks(idxs)
chunk_start_idxs = [];
chunk_end_idxs = [];

for i = 1:length(idxs)-1
    if i == 1
        chunk_start_idxs = [chunk_start_idxs, idxs(i)];
        continue;
    end

    if idxs(i+1) - idxs(i) ~= 1
        chunk_end_idxs = [chunk_end_idxs, idxs(i)];
        continue;
    end

    if idxs(i) - idxs(i-1) ~= 1
        chunk_start_idxs = [chunk_start_idxs, idxs(i)];
        continue;
    end
end

if length(chunk_start_idxs) ~= length(chunk_end_idxs)
   chunk_end_idxs = [chunk_end_idxs, idxs(end)]; 
end

end
