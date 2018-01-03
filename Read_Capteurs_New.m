close all
clear all
clc
clf

%base_folder = '/home/carvalhol/Desktop/Energy';
%base_folder = '/Users/carvalhol/Desktop/Energy';
%base_folder = '/Users/carvalhol/Desktop/Energy2';
base_folder = '/Users/carvalhol/Desktop/Energy4';
%base_folder = '/Users/carvalhol/Desktop/Energy_case_new';
%base_folder = '/Users/carvalhol/Desktop/res_Cineticas';

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

% 1      25.000000        25.000000        25.000000 
% 2      25.000000        25.000000        26.000000 
% 3      25.000000        25.000000        27.000000 
% 4      25.000000        25.000000        28.000000 
% 5      25.000000        26.000000        25.000000 
% 6      25.000000        27.000000        25.000000 
% 7      25.000000        28.000000        25.000000 
% 8      26.000000        25.000000        25.000000 
% 9      27.000000        25.000000        25.000000 
% 10     28.000000        25.000000        25.000000 
% 11     26.000000        26.000000        26.000000 
% 12     27.000000        27.000000        27.000000 
% 13     28.000000        28.000000        28.000000 
sensor = [2, 3, 4];
distance = [1, 2, 3];

%sensor(i) = 1; %Capteur in source

if (numel(capteurs) > 0)

    nTimeSteps = size(capteurs{1}.data, 1);
    nSensors = numel(sensor);
    
    u = zeros(nTimeSteps, 3, nSensors);
    a = zeros(nTimeSteps, 3, nSensors);
    t = zeros(nTimeSteps, 1, nSensors);
    
    for i = 1:numel(sensor)
         %Displacement
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Displ      1'));
         u(:,1,i) = capteurs{sensor(i)}.data(:,pos);
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Displ      2'));
         u(:,2,i) = capteurs{sensor(i)}.data(:,pos);
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Displ      3'));
         u(:,3,i) = capteurs{sensor(i)}.data(:,pos);
         %Acceleration
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Accel      1'));
         a(:,1,i) = capteurs{sensor(i)}.data(:,pos);
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Accel      2'));
         a(:,2,i) = capteurs{sensor(i)}.data(:,pos);
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Accel      3'));
         a(:,3,i) = capteurs{sensor(i)}.data(:,pos);
         %Time
         pos = find(strcmp(capteurs{sensor(i)}.labels, 'Time       1'));
         t(:,i) = capteurs{sensor(i)}.data(:,pos);
    end

%% Analitical Displacement
% Eduardo Kausel page 48
% Three-dimensional problems in full,homogeneous spaces
% 4.2Point load (Stokes problem)

alpha = 5; %Vp
beta = alpha/sqrt(2); %Vs
rho = 2800; %Density
mu = rho*beta^2; %G

T = 1.5; %period
nTimeSteps = 868;

%Sensor coordinates (relative to source)
R_mat=[1, 0, 0; ...
       2, 0, 0; ...
       3, 0, 0];

%Calculate displacement evolution for each R

u_signal=zeros(size(R_mat,1),nTimeSteps);

for coord=1:size(R_mat,1)

    R_vec = R_mat(coord,:);
    R = norm(R_vec);
    t_p=R/alpha;
    t_s=R/beta;


    %Time evolution
    u_Rx=zeros(1,nTimeSteps);
    count=0;
    t_vec=linspace(0,T,nTimeSteps);

    for count=1:nTimeSteps
        t_now=t_vec(count);
        dirac_p=1*(t_now==t_p);
        dirac_s=1*(t_now==t_s);
        hside_p=1/2*(t_now==t_p)+(t_now>t_p);
        hside_s=1/2*(t_now==t_s)+(t_now>t_s);

        PSI=dirac_s - t_now*(((beta/alpha)^2*hside_p/t_p^2)-hside_s/t_s^2);
        X=(beta/alpha)^2*dirac_p-dirac_s+3*t_now*(((beta/alpha)^2)*hside_p/t_p^2 -hside_s/t_s^2);

        u_Rx(count) = (PSI + X)/(4*pi*mu*R);

    end

    Signal = rickerSEM(t_vec,2,0.5,1);
    u_signal(coord,:) = ifft(fft(Signal).*fft(u_Rx));
end


    
%% Plot
figure (1)
%X
subplot(1,3,1);
hold on
count=0;
for i = 1:numel(sensor)
    count = count + 1;
    plot(t(:,i), u(:,1,i));
    label{i} = ['R = ',num2str(distance(i))];
end
for i = 1:size(u_signal,1)
    count = count + 1;
    plot(t_vec, u_signal(i,:));
    label{count} = ['Analytic R = ',num2str(distance(i))];
end
%title('Displacement X','FontSize',20);
legend(label,'FontSize',20);
xlabel('Time [s]','FontSize',20);
ylabel('Displacement [m]','FontSize',20);
set(gca,'fontsize',20)
hold off
%Y
subplot(1,3,2);
hold on
for i = 1:numel(sensor)
    plot(t(:,i), u(:,2,i));
    label{i} = ['R = ',num2str(distance(i))];
end
title('Displacement Y','FontSize',20);
legend(label,'FontSize',20);
xlabel('Time [s]','FontSize',20);
ylabel('Displacement [m]','FontSize',20);
set(gca,'fontsize',20)
hold off
%Z
subplot(1,3,3);
hold on
for i = 1:numel(sensor)
    plot(t(:,i), u(:,3,i));
    label{i} = ['R = ',num2str(distance(i))];
end
title('Displacement Z','FontSize',20);
legend(label,'FontSize',20);
xlabel('Time [s]','FontSize',20);
ylabel('Displacement [m]','FontSize',20);
set(gca,'fontsize',20)
hold off


%% 



    
%     %Energy
%     prop_size = 1;
%     En_P = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
%     pos = find(strcmp(capteurs{capt_nb}.labels, 'EnergyP    1'));
%     En_P(:) = capteurs{capt_nb}.data(:,pos);
%     prop_size = 1;
%     En_S = zeros(size(capteurs{capt_nb}.data, 1) , prop_size);
%     pos = find(strcmp(capteurs{capt_nb}.labels, 'EnergyS    1'));
%     En_S(:) = capteurs{capt_nb}.data(:,pos);
    
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

% %% Calculating Injected Energy
% 
% %Read case information
% %cd(base_folder);
% %c_folder = '/home/carvalhol/Desktop/Energy'; %Case folder, has material.input and input.spec
% %info = read_SEM_files(c_folder);
% %vs = info.Vs;
% %vp = info.Vp;
% %f = info.freq;
% %tau = info.tau;
% 
% %Making source
% t = E(:,1);
% f = 2;
% tau = 0.5;
% F = zeros(size(u));
% F(:,1) = -ricker(t, f, tau); %Ricker on X+
% 
% element_vol = (0.5)^3;
% dom_vol = (50)^3;
% 
% En_src=F.*u; %poid point de gauss coin
% En_src=F.*u*(1/10)^3*element_vol*8; %poid point de gauss coin
% 
% En_src_int = trapz(t,En_src);
% 
% %Cumulated Energy ?
% E_cum_src = zeros(numel(t)-1, 1);
% for i=1:numel(E_cum_src)
%     E_cum_src(i) = trapz(t(1:i+1),En_src(1:i+1));
% end
% 
% E_cum = zeros(numel(t)-1, 1);
% for i=1:numel(E_cum)
%     E_cum(i) = trapz(t(1:i+1),E(1:i+1,6));
% end
% 
% %Plot
% En_labels = {'P-Energy', 'S-Energy', 'Residual Energy', ...
%              'Kinetic Energy', 'Total Energy'};
% 
% %Source
% figure (1)
% subplot(2,2,1);
% title('Source (MATLAB)')
% hold on
% plot(t, F(:,1));
% plot(t, F(:,2));
% plot(t, F(:,3));
% legend({'Source x','Source y','Source z'},'FontSize',20)
% xlabel('Time [s]')
% ylabel('Force [N]')
% set(gca,'fontsize',20)
% hold off
% 
% %Displacement
% %figure (2)
% subplot(2,2,2);
% title('Displacement (Simulation)')
% hold on
% plot(t, u(:,1));
% plot(t, u(:,2));
% plot(t, u(:,3));
% legend({'Displacement x','Displacement y','Displacement z'},'FontSize',20)
% xlabel('Time [s]')
% ylabel('Displacement [m]')
% set(gca,'fontsize',20)
% hold off
% 
% %Energy
% %figure (3)
% subplot(2,2,3);
% title('Energy Source (WHAT WE WANT)')
% hold on
% plot(t, En_src(:,1));
% plot(t, En_src(:,2));
% plot(t, En_src(:,3));
% plot(t(1:end-1), E_cum_src(:))
% indexmax = find(max(E_cum_src) == E_cum_src);
% %xmax = t(indexmax);
% xmax = 4;
% ymax = E_cum_src(indexmax);
% strmax = ['Maximum = ',num2str(ymax)];
% text(xmax,ymax,strmax,'HorizontalAlignment','right','FontSize',20);
% legend({'Energy x','Energy y','Energy z','Cumulated Energy'},'FontSize',20)
% xlabel('Time [s]')
% ylabel('Energy [J]')
% set(gca,'fontsize',20)
% hold off
%          
% %figure (4)
% subplot(2,2,4);
% title('Energy Media (WHAT SIMULATION GIVES US)')
% hold on
% %plot(t, E(:,6))
% plot(t, E(:,pos_2)*element_vol);
% plot(t, E(:,pos_3)*element_vol);
% plot(t, E(:,pos_4)*element_vol);
% plot(t, E(:,pos_5)*element_vol);
% plot(t, E(:,pos_6)*element_vol);
% legend(En_labels,'FontSize',20)
% xlabel('Time [s]')
% ylabel('Energy [J]')
% set(gca,'fontsize',20)
% 
% %plot(t, E(:,6))
% %plot(t(1:end-1), E_cum(:))
% hold off
% 
% 
% % figure (4)
% % title('Energy Media ')
% % hold on
% % %plot(t, E(:,6))
% % plot(t, E(:,pos_2));
% % plot(t, E(:,pos_3));
% % plot(t, E(:,pos_4));
% % plot(t, E(:,pos_5)*10);
% % plot(t, E(:,pos_6)+9*E(:,pos_5));
% % legend(En_labels,'FontSize',20)
% % xlabel('Time [s]')
% % ylabel('Energy [J]')
% % set(gca,'fontsize',20)
% 
% %E_cum(end);



