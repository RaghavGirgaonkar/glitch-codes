function [] = createHDF5file(outfilename, whtndfiltdata, TF, N)
%Create HDF5 file for Whitened data
% Create HDF5 file
file_id = H5F.create(outfilename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

% h5create(outfilename,"/strain",[2 N]);
h5create(outfilename,"/strain/Strain",[N 1]);
h5create(outfilename,"/strain/condTF",[length(TF) 1]);
h5write(outfilename, '/strain/Strain', whtndfiltdata');
h5write(outfilename, '/strain/condTF', TF');

% Set attributes
strainDS = h5info(outfilename, '/strain/Strain');
strainDS.Attributes.Name = 'Whitened Strain';

condTFDS = h5info(outfilename, '/strain/condTF');
condTFDS.Attributes.Name = 'Transfer Function';

% Close file
H5F.close(file_id);


end

