function [info] = read_SEM_files(caseFolder)
    
    root = pwd;
    cd(caseFolder);
    info.caseFolder = caseFolder;
    
    fileText = fileread('input.spec');
    
    info.freq = findInfo('freq =', '%f', fileText);
    info.func = findInfo('func =', '%s', fileText);
    info.tau  = findInfo('tau =', '%s', fileText);
    info.traces_format = findInfo('traces_format =', '%s', fileText);
    info.capteurs_file = 'stations.txt';
    info.outVar_Text = {'enP', 'enS', 'evol', 'pre', 'dis', ...
                        'vel', 'acc', 'edev', 'sdev'};
    outVar_Dim = int32([1, 1, 1, 1, 3, 3, 3, 6, 6]');
    info.outVar_Int  = zeros(numel(info.outVar_Text), 1, 'int32');
    
    for i =1:numel(info.outVar_Text)
        info.outVar_Int(i) = findInfo([info.outVar_Text{i},' ='], '%d', fileText);
    end
    
    info.outVar_Int = info.outVar_Int .* outVar_Dim;
    
    fid = fopen('material.input', 'r' );
    nDomains = fscanf(fid,'%d');
    for i = 1:nDomains
        data = fscanf(fid,'%s %g %g %g %d %g %g',[1 7]);
        mat_type = char(data(1));
        if(strcmp(mat_type, 'S') || strcmp(mat_type, 'R'))
            info.Vp = data(2);
            info.Vs = data(3);
            info.mat_nb = i;
            info.mat_type = mat_type;
            break;
        end
    end
    
    fclose(fid);
    
    cd(root);

end

function [founded_info] = findInfo(pattern, type, fileText)
    
    %text = fileread(file);
    index = strfind(fileText, pattern);
    beg = index(1) + numel(pattern);
    temp = textscan(fileText(beg:end),type);
    founded_info = temp{1};
    if(strcmp(type,'%s'))
        founded_info = strjoin(founded_info(1));
        
        if(strcmp(founded_info(end),';'))
            %founded_info = textscan(founded_info(1:end-1),type);
            %founded_info = strjoin(founded_info(1));
            founded_info = founded_info(1:end-1);
        end
    end
    
end