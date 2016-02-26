clc
clear
close all

%% USER
%baseFolder = '/home/lcp/Desktop/RF_Matlab';
%baseFolder = '/mssmat2/home/paludo/Desktop/RF_Matlab';
baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';
searchFolder = '2DSamples/2D/FFT-g';
fileName = 'BBox_L01_001_001.h5';
nFolders = 1;

for i = 1:nFolders
    searchFolderIt = strcat(searchFolder,'/', sprintf('%03d',i));
    cd([searchFolder])
    pathList = subdir(fileName); %Gives all the paths that have a singleGen file
    cd(baseFolder);
    vecSize = size(pathList,1);
    
    Sk_1  = cell(vecSize,1);
    Sk_2  = cell(vecSize,1);
    
    for file = 1:vecSize
        info = hdf5info(pathList(file).name);
        nAttrib = size(info.GroupHierarchy.Attributes,2);
        nDSet  = size(info.GroupHierarchy.Datasets,2);
        
        for j = 1 : nAttrib
            if (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Sk_1'))
                Sk_1{file} = info.GroupHierarchy.Attributes(j).Value;
            elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Sk_2'))
                Sk_2{file} = info.GroupHierarchy.Attributes(j).Value;
            end
        end
        
    end
end

