close all
clear all

c_folder = './Energy_Results'; %Case folder, has material.input and input.spec
t_folder = './Energy_Results'; %Traces Folder

%% Read traces
cd(t_folder)
files = dir('capteurs*');
capteurs = cell(0);

for i = 1 : numel(files)
    
    info = hdf5info(files(i).name);
    labels = h5read(files(i).name, '/Variables');
    
    for j = 1 : size(info.GroupHierarchy.Datasets,2)
        
        if(strcmp(info.GroupHierarchy.Datasets(j).Name, '/En_PS') || ...
           strcmp(info.GroupHierarchy.Datasets(j).Name, '/Variables'))
            continue;
        end
        
        data = h5read(files(i).name, info.GroupHierarchy.Datasets(j).Name);
        data = data(1:numel(labels), :)';
        
        capteurs{end+1} = struct('labels', {labels}, ...
                                  'data', data, ...
                                  'dset', info.GroupHierarchy.Datasets(j).Name, ...
                                  'file', files(i).name); 
       
        clear data
    end
end

%% Read case information
info = read_SEM_files(c_folder);
f = info.freq;
tau = info.tau;
vs = info.Vs;
vp = info.Vp;

%% Take the information you want



