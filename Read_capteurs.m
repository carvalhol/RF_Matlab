close all
clear all
clc
clf


base_folder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';

%% Read traces
cd(base_folder);
t_folder = './Energy_Results'; %Traces Folder
cd(t_folder)
files = dir('capteurs*');
capteurs = cell(0);
capteurs_En = cell(0);

for i = 1 : numel(files)
    
    info = hdf5info(files(i).name);
    labels = h5read(files(i).name, '/Variables');
    
    for j = 1 : size(info.GroupHierarchy.Datasets,2)
        
        if(strcmp(info.GroupHierarchy.Datasets(j).Name, '/Variables'))
            continue;
        end
        
        data = h5read(files(i).name, info.GroupHierarchy.Datasets(j).Name);
                
        if(strcmp(info.GroupHierarchy.Datasets(j).Name, '/En_PS'))
            data = data';
            capteurs_En{end+1} = struct('labels', {{'Time        '; ...
                                                'En_P        '; ...
                                                'En_S        '; ...
                                                'En_R        '; ...
                                                'En_Cin      '; ...
                                                'En_Tot      ' ...
                                                }}, ...
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

%Interpolation Sensors
capt_nb = 1;

if (numel(capteurs) > 0)
    %Displacement
    prop_size = 3;
    u = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
    pos = find(strcmp(capteurs{capt_nb}.labels, 'Displ      1'));
    u(:,1) = capteurs{capt_nb}.data(:,pos);
    pos = find(strcmp(capteurs{capt_nb}.labels, 'Displ      2'));
    u(:,2) = capteurs{capt_nb}.data(:,pos);
    pos = find(strcmp(capteurs{capt_nb}.labels, 'Displ      3'));
    u(:,3) = capteurs{capt_nb}.data(:,pos);

    %Time
    prop_size = 1;
    t = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
    pos = find(strcmp(capteurs{capt_nb}.labels, 'Time        '));
    t(:) = capteurs{capt_nb}.data(:,pos);
end

%Energy Sensors
if (numel(capteurs_En) > 0)
    %Enegies
    prop_size = 5;
    E = zeros(size(capteurs_En{1}.data, 1) , prop_size);
    pos = find(strcmp(capteurs_En{1}.labels, 'En_P        '));
    E(:,1) = capteurs_En{1}.data(:,pos);
    pos = find(strcmp(capteurs_En{1}.labels, 'En_S        '));
    E(:,2) = capteurs_En{1}.data(:,pos);
    pos = find(strcmp(capteurs_En{1}.labels, 'En_R        '));
    E(:,3) = capteurs_En{1}.data(:,pos);
    pos = find(strcmp(capteurs_En{1}.labels, 'En_Cin      '));
    E(:,4) = capteurs_En{1}.data(:,pos);
    pos = find(strcmp(capteurs_En{1}.labels, 'En_Tot      '));
    E(:,5) = capteurs_En{1}.data(:,pos);
    %Time
    prop_size = 1;
    t_En = zeros(size(capteurs_En{capt_nb}.data, 1) , prop_size);
    pos = find(strcmp(capteurs_En{capt_nb}.labels, 'Time        '));
    t_En(:) = capteurs_En{capt_nb}.data(:,pos);

end

%Plot
figure (1)
plot(t, u(:,1))


%% Calculating Injected Energy

%Read case information
%cd(base_folder);
%c_folder = './Energy_Results'; %Case folder, has material.input and input.spec
%info = read_SEM_files(c_folder);
%vs = info.Vs;
%vp = info.Vp;
%f = info.freq;
%tau = info.tau;

%Making source
%f = 3;
%tau = 0.4;
%source = ricker(t, f, tau);

%Plot
%figure (2)
%plot(t, source(:,1))

du = circshift(u(:,1),1)-u(:,1);
%F  = (circshift(source(:,1),1)+source(:,1))/2;
delta_t = t(2)-t(1);
du_dt = du(1:end-1)/delta_t;
%energy = trapz(t(1:end-1),source(1:end-1,1).*du_dt(:,1))
%energy_abs = trapz(t(1:end-1),abs(source(1:end-1,1).*du_dt(:,1)));

energy_Fdu = abs(F(1:end-1)'*du(1:end-1));


%energy = trapz(t,abs(source(:,1).*u(:,1)))


%Cumulated Energy
E_cum = zeros(numel(t_En)-1, 1);
for i=1:numel(t_En)-1
    E_cum(i) = trapz(t_En(1:i+1),E(1:i+1,5));
end

%Plot
figure (3)
hold on
plot(t_En, E(:,5))
plot(t_En(1:end-1), E_cum(:))
hold off

E_cum(end)



