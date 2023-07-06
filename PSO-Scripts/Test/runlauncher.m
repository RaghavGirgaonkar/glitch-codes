hdf5filenames   = 'hdf5_file_names.txt';
segfile = 'segments.txt';
psdfile = 'psds.mat';
sampFreq = 4096;
for i = 1:10
    launcherscript(i, hdf5filenames,segfile,psdfile, sampFreq);
end