close all
clear all

%cd ./GLOBAL_FIELD/h5
%cd ./GLOBAL_PROC/h5
%cd ./INDEP_FIELD/h5
%cd ./INDEP_PROC/h5
cd /Users/carvalhol/Desktop

plotRes = true;
init = pwd;

nome = dir('*.h5');

for i = 1 : size(nome,1)
    info = hdf5info(nome(i).name);
    nAttrib = size(info.GroupHierarchy.Attributes,2);
    nDSet  = size(info.GroupHierarchy.Datasets,2);
    
    %Reading Attributes
    for j = 1 : nAttrib
        if (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'independent'))
            independent = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'nb_procs'))
            nb_procs = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'nDim'))
            nDim = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Nmc'))
            Nmc = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'method'))
            method = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seedStart'))
            seedStart = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'corrMod'))
            corrMod = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'margiFirst'))
            margiFirst = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_CPU_Time'))
            gen_CPU_Time = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_WALL_Time'))
            gen_WALL_Time = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seed'))
            seed = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xMaxGlob'))
            xMaxGlob = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xMinGlob'))
            xMinGlob = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xStep'))
            xStep = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'corrL'))
            corrL = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'overlap'))
            overlap = info.GroupHierarchy.Attributes(j).Value;
        end
    end
    
    
    for j = 1 : nDSet
        RF = h5read(nome(i).name, info.GroupHierarchy.Datasets(j).Name);
    end
    
    xNStep = (xMaxGlob - xMinGlob)./xStep + 1;
    RF = reshape(RF, xNStep');
    
    
    %     for j = 1 : size(info.GroupHierarchy.Datasets,2)
    %         aux{i,j} = h5read(nome(i).name, info.GroupHierarchy.Datasets(j).Name);
    %     end
    
%     %% Ploting Results
%     if(plotRes)
%         switch nDim
%             case 1
%                 
%             case 2
%                 %figure(1)
%                 %hold on
%                 [x,y] = meshgrid(xMinGlob(1):xStep(1):xMaxGlob(1),xMinGlob(2):xStep(2):xMaxGlob(2));
%                 maxPlot = 1;
%                 hSurf = surf(x,y,RF,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%                 %hold off
%             case 3
%                 
%             otherwise
%                 error('Dimension not accepted in plot')
%         end
%     end
end

%cd init
cd ../..


%% Statistics
dir = double(1:1:nDim);
C = cell(nDim,1);
lc = zeros(nDim,1);

for i = 1: nDim
    [C{i}, lc(i)] = corr_RF(RF, RF, dir(i), xStep(i));
    %lc(i) = 2*trapz(C{i})*xStep(i)*(2*numel(C{i})-1)/(2*numel(C{i}));
end

avg    = mean(RF(:));
stdDev = std(RF(:));

%clf
%plot(xStep(1)*(1:numel(C{1}))*(numel(C{1})-1)/numel(C{1}),C{1});
%Ctot = (C{1}+C{2})/2;
%lcTot = 2*trapz(Ctot)*xStep(1)*(2*numel(Ctot)-1)/(numel(Ctot)*2);
%plot(xStep(1)*(1:numel(Ctot)),Ctot);