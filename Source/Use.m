delete('warnings.log');
cd('../Data');
filenames = dir('*.mat');
for kk = 1:numel(filenames)
    disp(filenames(kk).name);
    load(filenames(kk).name);
    cd('../Source');
    lng=dataTable_gps(:,2);
   lat=dataTable_gps(:,3);    
    [eas,nor,ele] = OSTN15_Matlab('gps-to-grid',lat,lng);
   cd('../SortedData');
   savenew = strcat(filenames(kk).name(:,1:8), "_osgb-elev", filenames(kk).name(:,13:end-4), ".mat");
  save(savenew, 'eas', 'nor', 'ele');
   cd('../Data');
   clearvars('-except', 'kk', 'filenames');
end