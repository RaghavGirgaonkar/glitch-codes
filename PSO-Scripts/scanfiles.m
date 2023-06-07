function []=scanfiles(dirName)

%Make List of all HDF5 files in specified directory
fileList = dir(fullfile(dirName, '*.hdf5'));

%Print Filenames
for i = 1:numel(fileList)
    fprintf('%s\n', fileList(i).name);
end

%Get GPS Start times and duration of each file and store in separate txt
%file

%Open txt file 
outputFile = fullfile(dirName, 'hdf5_file_names.txt');
fid = fopen(outputFile, 'w');

for i = 1:numel(fileList)
    name = fileList(i).name;
    GPSstart = h5read(name,'/meta/GPSstart');
    duration = h5read(name, '/meta/Duration');
    fprintf(fid, '%s\t%d\t%d\n', name,GPSstart,GPSstart+duration);
end

fclose(fid);

end