close all
clear all
clc
clf

base_folder = '/home/carvalhol/Desktop/Energy';

%% Read traces
cd(base_folder);
%t_folder = './Energy_Results'; %Traces Folder
%t_folder = '/home/carvalhol/Desktop/traces';
t_folder = 'traces';
cd(t_folder)
files = dir('capteurs*.h5');
capteurs = cell(0);
capteurs_En = cell(0);

%% FINDING LABELS
for i = numel(files) : -1: 1 
    
    files(i).name
    
    if(~endsWith(files(i).name,'.h5'))
        files(i).name = {};
    end
    
    info = hdf5info(files(i).name);
    
    for j = 1 : size(info.GroupHierarchy.Datasets,2)
        if(strcmp(info.GroupHierarchy.Datasets(j).Name, '/Variables'))
            labels = h5read(files(i).name, '/Variables');
            continue;
        elseif(strcmp(info.GroupHierarchy.Datasets(j).Name, '/En_PS_Variables'))
            Eng_labels = h5read(files(i).name, '/En_PS_Variables');
            continue;
        end
    end
    
end

for i = 1 : numel(files)
    
    info = hdf5info(files(i).name);
    
    for j = 1 : size(info.GroupHierarchy.Datasets,2)
        
        if(endsWith(info.GroupHierarchy.Datasets(j).Name,'_pos') || ...
           endsWith(info.GroupHierarchy.Datasets(j).Name,'_Variables'))
            continue;
        end      
        
        data = h5read(files(i).name, info.GroupHierarchy.Datasets(j).Name);
                
        if(strcmp(info.GroupHierarchy.Datasets(j).Name, '/En_PS'))
            data = data';
            capteurs_En{end+1} = struct('labels', {Eng_labels}, ...
                                      'data', data, ...
                                      'dset', info.GroupHierarchy.Datasets(j).Name, ...
                                      'file', files(i).name); 
        else
            data = data(1:numel(labels), :)';

            capteurs{end+1} = struct('labels', {labels}, ...
                                      'data', data, ...
                                      'dset', info.GroupHierarchy.Datasets(j).Name, ...
                                      'file', files(i).name);
        end
        clear data
    end
end

%% Take the information you want

source_capt = 1; %Capteur in source
if (numel(capteurs) > 0)
     %Displacement
     prop_size = 3;
     u = zeros(size(capteurs{source_capt}.data, 1) , prop_size);
     a = zeros(size(capteurs{source_capt}.data, 1) , prop_size);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Displ      1'));
     u(:,1) = capteurs{source_capt}.data(:,pos);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Displ      2'));
     u(:,2) = capteurs{source_capt}.data(:,pos);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Displ      3'));
     u(:,3) = capteurs{source_capt}.data(:,pos);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Accel      1'));
     a(:,1) = capteurs{source_capt}.data(:,pos);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Accel      2'));
     a(:,2) = capteurs{source_capt}.data(:,pos);
     pos = find(strcmp(capteurs{source_capt}.labels, 'Accel      3'));
     a(:,3) = capteurs{source_capt}.data(:,pos);
 
%     %Time
%     prop_size = 1;
%     t = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
%     pos = find(strcmp(capteurs{capt_nb}.labels, 'Time       1'));
%     t(:) = capteurs{capt_nb}.data(:,pos);
%     
% %     %Energy
% %     prop_size = 1;
% %     En_P = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
% %     pos = find(strcmp(capteurs{capt_nb}.labels, 'EnergyP    1'));
% %     En_P(:) = capteurs{capt_nb}.data(:,pos);
% %     prop_size = 1;
% %     En_S = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
% %     pos = find(strcmp(capteurs{capt_nb}.labels, 'EnergyS    1'));
% %     En_S(:) = capteurs{capt_nb}.data(:,pos);
%     
end

%Energy Sensors
prop_size = numel(capteurs_En{1}.labels);
E = zeros(size(capteurs_En{1}.data, 1) , prop_size);
pos_1 = find(strcmp(capteurs_En{1}.labels, 'Time       1'));
pos_2 = find(strcmp(capteurs_En{1}.labels, 'Eng_P      1'));
pos_3 = find(strcmp(capteurs_En{1}.labels, 'Eng_S      1'));
pos_4 = find(strcmp(capteurs_En{1}.labels, 'Eng_Resid  1'));
pos_5 = find(strcmp(capteurs_En{1}.labels, 'Eng_Cine   1'));
pos_6 = find(strcmp(capteurs_En{1}.labels, 'Eng_Total  1'));

E(:,1) = capteurs_En{i}.data(:,pos_1);

for i = 1 : numel(capteurs_En)
    E(:,2) = E(:,2) + capteurs_En{i}.data(:,pos_2);   
    E(:,3) = E(:,3) + capteurs_En{i}.data(:,pos_3);    
    E(:,4) = E(:,4) + capteurs_En{i}.data(:,pos_4);    
    E(:,5) = E(:,5) + capteurs_En{i}.data(:,pos_5);    
    E(:,6) = E(:,6) + capteurs_En{i}.data(:,pos_6);
end


%Plot
%figure (1)
%plot(t, u(:,1))
%plot(E(:,1), E(:,6))

%% Calculating Injected Energy

%Read case information
%cd(base_folder);
%c_folder = '/home/carvalhol/Desktop/Energy'; %Case folder, has material.input and input.spec
%info = read_SEM_files(c_folder);
%vs = info.Vs;
%vp = info.Vp;
%f = info.freq;
%tau = info.tau;

%Making source
t = E(:,1);
f = 2;
tau = 0.5;
F = zeros(size(u));
F(:,1) = ricker(t, f, tau); %Ricker on X+
A = ricker(t, f, tau);

En_src=F.*u;

En_src_int = trapz(t,En_src);


%Cumulated Energy ?
E_cum = zeros(numel(t)-1, 1);
for i=1:numel(E_cum)
    E_cum(i) = trapz(t(1:i+1),E(1:i+1,6));
end

%Plot
En_labels = {'P-Energy', 'S-Energy', 'Residual Energy', ...
             'Cinetic Energy', 'Total Energy'};
figure (3)
hold on
%plot(t, E(:,6))
plot(t, E(:,pos_2));
plot(t, E(:,pos_3));
plot(t, E(:,pos_4));
plot(t, E(:,pos_5));
plot(t, E(:,pos_6));
legend(En_labels)

%plot(t, E(:,6))
%plot(t(1:end-1), E_cum(:))
hold off

%E_cum(end);



