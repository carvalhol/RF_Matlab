clc
clear
close all

cleaning=false;
launch_folder='/Users/carvalhol/Desktop/WEAK';
%launch_folder='/home/carvalhol/Desktop/WEAK';
cd(launch_folder);

%%CONSTANTS
fileName = 'INFO-Sample_1.h5';
f_res_name = 'WEAK_data.mat';

%% USER
orig_folder = pwd;
paths = {...
    '/Users/carvalhol/Desktop/WEAK/3D/FFT-l'...
}; %Path till iterations subfolders
%paths = {...
%    '/home/carvalhol/Desktop/WEAK/3D/FFT-l'...
%}; %Path till iterations subfolders



%% PROCESSING

nCases = numel(paths);
res = cell(nCases,1);

f_list_name = '_list.txt';
f_list_name = [orig_folder,'/',f_list_name];

if(exist(f_list_name,'file')==2 && cleaning)
    [status, probFile] = system(['tail -n 1 ',f_list_name]);
    probFile = strtrim(probFile);
    %probFile = string(probFile);
    disp(['Last file on list: ', probFile])
    prompt = 'Delete last folder (Y/n) [n]: ';
    ans = input(prompt,'s');
    if isempty(ans)
        ans = 'n';
    end
    if(strcmp(ans,'Y'))
        if(exist(probFile,'file')==2)
            delete(probFile)
        else
           disp('File does not exist') 
        end
    end
    delete(f_list_name)
end

f_res_path = [orig_folder,'/',f_res_name];
if~(exist(f_res_path,'file')==2)

    fileID = fopen(f_list_name,'w');

    for p = 1:nCases

        cd(orig_folder);
        cd(paths{p});
        foundFiles = subdir(fileName);
        nIter=numel(foundFiles);

        %Things to read just once inside the path
        f_name = foundFiles(1).name;
        res{p}.time_labels = h5read(f_name, '/time_labels');

        res{p}.iter = cell(nIter,1);

        for it = 1:nIter
            [folder,file,ext] = fileparts(foundFiles(it).name);
            cd(folder);
            files = dir('INFO-*.h5');
            nSamples=numel(files);

            %Read just once inside each iteration
            f_name = files(1).name;
            res{p}.iter{it}.nb_procs = h5readatt(f_name,'/','nb_procs');
            res{p}.iter{it}.corrL = h5readatt(f_name,'/','corrL');
            res{p}.iter{it}.xMaxGlob = h5readatt(f_name,'/','xMaxGlob');
            res{p}.iter{it}.xMinGlob = h5readatt(f_name,'/','xMinGlob');
            res{p}.iter{it}.overlap = h5readatt(f_name,'/','overlap');
            res{p}.iter{it}.nFields = h5readatt(f_name,'/','nFields');
            res{p}.iter{it}.L = h5readatt(f_name,'/','L');
            res{p}.iter{it}.Np = h5readatt(f_name,'/','Np');
            res{p}.iter{it}.Np_ovlp = h5readatt(f_name,'/','Np_ovlp');

            times = zeros(numel(res{p}.time_labels),...
                          res{p}.iter{it}.nb_procs,...
                          numel(files));

            for f = 1:nSamples

                f_name = files(f).name;
                fprintf(fileID,[folder,'/',f_name,'\n']);
                times(:,:,f) = h5read(f_name, '/times');
            end

            res{p}.iter{it}.times_max = max(max(times,[],3),[],2);
            res{p}.iter{it}.times_min = min(min(times,[],3),[],2);
            res{p}.iter{it}.times_avg = mean(mean(times,3),2);
            res{p}.iter{it}.times_std = mean(std(times,1,2),3);
            %times_tmp = times(9,:,:); %Wall time
            times_tmp = sum(times(1:8,:,:),1);
            times_tmp = sum(times_tmp,2);
            res{p}.iter{it}.times_plot = times_tmp(:);

        end

    end
    save(f_res_path,'res')
    fclose(fileID);
else
    load(f_res_path)
end

%% Plot

figure(1)
for p = 1:nCases
    
    nIter = numel(res{p}.iter);
    x  = zeros(nIter, 1);
    %y1 = zeros(nIter, 1);
    %y2 = zeros(nIter, 1);
    %y3 = zeros(nIter, 1);
    
    grp = [];
    y=[];
    labels = cell(1,nIter);
    
    for it = 1:nIter
        dx_lc = (res{p}.iter{it}.xMaxGlob - res{p}.iter{it}.xMinGlob)...
                 ./res{p}.iter{it}.corrL;
        x(it) = dx_lc(1);
        %y1(it) = sum(res{p}.iter{it}.times_avg(1:8));
        y = [y,res{p}.iter{it}.times_plot'];
        grp = [grp, x(it)*ones(1,numel(res{p}.iter{it}.times_plot))];
        nProcs = res{p}.iter{it}.nb_procs;
        %boxplot(y,x)
        %y2(it) = sum(res{p}.iter{it}.times_min(1:7));
        %y3(it) = sum(res{p}.iter{it}.times_max(1:7));
        %y2(it) = y1(it) + sum(res{p}.iter{it}.times_std(1:8));
        %y3(it) = y1(it) - sum(res{p}.iter{it}.times_std(1:8));
        
        labels(it) = [strcat(num2str(x(it)),{'lc '})];
        %labels(it) = [strcat(num2str(x(it)),{'lc '}); ...
        %              strcat(num2str(nProcs),{'procs'})];
        %labels{it} = sprintf('%d lc \n %d procs', x(it),nProcs);
        %labels{it} = '\begin{tabular}{c} line 1 \\ line 2 \\ line 3\\ line 4\end{tabular}';
    end

    y1 = mean(res{p}.iter{1}.times_plot);
    x1 = x(1);
    lin_scaling = y1 + y1*((x-x1)/x1).^3;
end

%hand = boxplot(y,grp,'Notch','on','Widths',0.5,...
%               'Labels',labels);
A={'1 2','2','3','4','5','6','7','8'};
hand = boxplot(y,grp,'Notch','on','Widths',0.5,'Labels',labels);
%set(gca,'xtick', 5, 'XTickLabel', labels, 'TickLabelInterpreter', 'latex')
set(hand(7,:),'Visible','off') 
hold on
plot(1:nIter,lin_scaling')
hold off

grid('on')
box('on')
%set(gca,'xscale','log')
set(gca,'yscale','log', 'FontSize',15)
% % legend(legendInfo(legendMask),'Location','northwest','FontSize',15)

%plot(x,lin_scaling);
%hold on
%hold off

%% TOCO

%data = rand(100,3);
%data = data(:);
%group = [1*ones(1,100),2*ones(1,100),3*ones(1,100)]';
%boxplot(data,group)
%boxplot(y,grp)
%hold on
%plot(1:3,[1,1,1])
%plot(1:nIter,lin_scaling')
%line(1:3,[1,1,1])

% % grid('on')
% % box('on')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log', 'FontSize',15)
% % legend(legendInfo(legendMask),'Location','northwest','FontSize',15)
% % % 
% % % saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');

% plot(x,y1)
% 
% plot_MaxMin = 1;
% 
% if(plot_MaxMin == 1)
%     grey = [0.4,0.4,0.4];
%     polygonX = [x; x(end:-1:1)];
%     polygonY = [y2; y3(end:-1:1)];
%     polygonZ = 0.0*polygonX -0.1;
%     %polygon = fill(polygonX, polygonY, grey);
%     polygon =  patch(polygonX, polygonY, polygonZ, grey);
%     set(polygon,'facealpha',.35);
%     set(get(get(polygon,'Annotation'),'LegendInformation'),...
%            'IconDisplayStyle','off'); % Exclude line from legend
% end


% %% CONSTANTS AND INITIALIZATION
% 
% %Constants
% nMethods = 4;
% kAdjust = 1;
% fileName = 'INFO-Sample_1.h5';
% resultsFolder = '/SAMPLES/';
% sizeBTVec = 8;
% nTests = numel(dims)*numel(methodStr);
% Legend=cell(0,1);
% plotCount = 0;
% set(0,'defaulttextinterpreter','latex')
% 
% linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
% '-.','--','-',':','-.','--','-',':','-.'));

%% LISTING FILES

% cd(baseFolder);
% cd(searchFolder);
% 
% folderList = subdir(fileName);
% vecSize = numel(folderList);
% 
% info = hdf5info(folderList(1).name);
% nAttrib = size(info.GroupHierarchy.Attributes,2);
% nDSet  = size(info.GroupHierarchy.Datasets,2);
% 
% 
% %% READING FILES AND PLOTTING DATA
% fig = figure(1);
% hold on
% hold all
% 
% for Idim = 1:numel(dims)
%     for Imet = 1:numel(methodStr)
%         
%         %Allocation---------------------
%         nb_procs = zeros(vecSize,1);
%         nDim = zeros(vecSize,1);
%         method = zeros(vecSize,1);
%         corrMod = zeros(vecSize,1);
%         gen_WALL_Time = zeros(vecSize,1);
%         nFields  = zeros(vecSize,dims(Idim));
%         procExtent = zeros(vecSize,dims(Idim));
%         xMaxGlob = zeros(vecSize,dims(Idim));
%         xMinGlob = zeros(vecSize,dims(Idim));
%         xStep    = zeros(vecSize,dims(Idim));
%         BT_min  = zeros(vecSize,sizeBTVec);
%         BT_max  = zeros(vecSize,sizeBTVec);
%         BT_avg  = zeros(vecSize,sizeBTVec);
%         BT_stdDev  = zeros(vecSize,sizeBTVec);
%         corrL   = zeros(vecSize,dims(Idim));
%         overlap = zeros(vecSize,dims(Idim));
%         localizationLevel = zeros(vecSize,1);
% 
%         for Ipath = 1:numel(folderList)
%             
%             %disp(folderList(Ipath).name)
%             info = hdf5info(folderList(Ipath).name);
%             nAttrib = size(info.GroupHierarchy.Attributes,2);
%             nDSet  = size(info.GroupHierarchy.Datasets,2);
%             %file_id = H5F.open(filename);
%             %H5F.close(file_id);
%             
%             %% Reading Attributes
%             for Iattr = 1 : nAttrib
%                 if (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nb_procs'))
%                     nb_procs(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nDim'))
%                     nDim(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
%                     %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Nmc'))
%                     %             Nmc(i) = info.GroupHierarchy.Attributes(j).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'method'))
%                     method(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
%                     %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seedStart'))
%                     %             seedStart(i) = info.GroupHierarchy.Attributes(j).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'corrMod'))
%                     corrMod(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
%                     %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'margiFirst'))
%                     %             margiFirst(i) = info.GroupHierarchy.Attributes(j).Value;
%                     %         %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_CPU_Time'))
%                     %         %    gen_CPU_Time(i) = info.GroupHierarchy.Attributes(j).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'gen_WALL_Time'))
%                     gen_WALL_Time(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nFields'))
%                     nFields(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'procExtent'))
%                     procExtent(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                     %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seed'))
%                     %    seed = info.GroupHierarchy.Attributes(j).Value;
%                  elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xMaxGlob'))
%                      xMaxGlob(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                  elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xMinGlob'))
%                      xMinGlob(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xStep'))
%                     xStep(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_min'))
%                     BT_min(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_max'))
%                     BT_max(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_avg'))
%                     BT_avg(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_stdDev'))
%                     BT_stdDev(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'corrL'))
%                     corrL(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'overlap'))
%                     overlap(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
%                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'localizationLevel'))
%                     localizationLevel(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;    
%                 end
%             end
%         end
%         
%         %% Cleaning vectors from repeated and empty files cases
%         
%         testMask  = false(vecSize, 1);
%         VT = false(vecSize, 1); %valid Tests
%         
%         labels = (1:vecSize)';
%         
%         if(time_Type == Wall_Time)
%             ylabel('Wall Time [s]', 'FontSize', 20);
%             time = sum(BT_avg(:,1:4),2); %Only Generation Time
%             %time = sum(BT_avg(:,1:8),2); %Whole Time (includi file writing)
%             %time = gen_WALL_Time;
%         elseif (time_Type == CPU_Time)
%             ylabel('CPU Time [s]', 'FontSize', 20);
%             time = sum(BT_avg(:,1:4),2).*nb_procs; %Only Generation Time
%             %time = sum(BT_avg(:,1:8),2).*nb_procs; %Whole Time (includi file writing)
%             %time = gen_WALL_Time;
%         end
%         
%         
%         %labels(time <= 0.0) = -1;
%         
%         for file = 1:vecSize
% 
%             %Testing if this test is not empty
%             if(time(file) > 0.0)
%                 VT(file) = true;
%             end
%             %Testing if this test is not redundant
%             if(labels(file) ~= file)
%                 continue
%             else
%                 testMask(file) = true;
%             end
%             
%             n = strfind(folderList(file).name, resultsFolder) - 1;
%             
%             str1 = folderList(file).name;
%             
%             for file2 = 1:vecSize
%                 if(file == file2)
%                     continue
%                 end
%                 
%                 str2 = folderList(file2).name;
%                 
%                 if(strncmpi(str1,str2,n))
%                     %'They Are In The Same Folder, average them'
%                     labels(file2) = file;
%                 end
%             end
%         end
%         
%         %Taking the time average and variance
%         var_time = accumarray(labels(VT), time(VT), [], @(x) var(x,1));
%         time_max = accumarray(labels(VT), time(VT), [], @max);
%         time_min = accumarray(labels(VT), time(VT), [], @min);
%         time_avg = accumarray(labels(VT), time(VT), [], @(x) mean(x,1));
%         toto = 4;
%         
%         %Number of points per proc
%         nPointsPerProc = zeros(numel(time_avg), 1);
%         pointsPerCorrL = 5;
%         extent = procExtent(testMask,:);
%         if(strcmp(methodStr{Imet}(end), 'g'))
%             extent = (xMaxGlob(testMask,:)-xMinGlob(testMask,:)); 
%             extent(:,Idim) = extent(:,Idim)./nb_procs(testMask);
%         end
%         pointsPerProc = prod(pointsPerCorrL*(extent./corrL(testMask,:)),2);
%         
%         disp(methodStr{Imet});
%         disp(dims(Idim));
%         disp('nPoints per Proc');
%         disp(pointsPerProc);
%         disp('Size Total (X, Y, Z)');
%         disp(xMaxGlob(testMask,:)-xMinGlob(testMask,:));
%       
%         
% 
%         %% Plotting Results
%         
%         lWidth = 2; %Line Width
%         mArea = 100; %MarkerSize
%         mSize = 9;
%         if(strcmp(methodStr{Imet}(end), 'l'))
%             plotId = 1;
%         else
%             plotId = 2;
%         end
%         plotId2 = dims(Idim);
%         makerStyle = {'o'; '^'; 's'};
%         makerFill = {'k'; 'k'; 'g'};
%         makerLine = {'k'; 'k'; 'g'};
%         colorStyle = {'k'; 'k'; 'g'};
%         lineStyle = {'o-'; '^-'};
%         if(plotId2 == 3)
%             lineStyle = {'o-.'; '^-.'};
%         end
%         
%         thisMarker = makerStyle{mod(plotId-1, size(makerStyle,1))+1};
%         thisFill   = makerFill{mod(plotId-1, size(makerFill,1))+1};
%         thisMLine  = makerLine{mod(plotId-1, size(makerLine,1))+1};
%         thisColor  = colorStyle{mod(plotId-1, size(colorStyle,1))+1};
%         thisLine   = lineStyle{mod(plotId-1, size(lineStyle,1))+1};
%                      
%         %Filtering
%         nPoints = prod(1+int64((xMaxGlob(testMask,:) - xMinGlob(testMask,:))./xStep(testMask,:)),2);
%         xRange = prod((xMaxGlob(testMask,:) - xMinGlob(testMask,:))./corrL(testMask,:),2);
%         nb_procs = nb_procs(testMask);
%         
%         
%         %Separe Number of processors
%         trial_processors = 1;
%         p_label = 1;
%         p_id    = ones(numel(nb_procs),1);
%         
%         if (trial_processors)
%             p_label = unique(nb_procs);
%             for i = 1:numel(p_label)
%                 p_id(nb_procs==p_label(i))=i;
%             end
%         end
%         
%         for i =1:numel(p_label)
%             plotCount = plotCount+1;
%             procMask = testMask(testMask);
%             procMask = (procMask & nb_procs==p_label(i));
%             
%             Legend(end+1) = strcat(testType,{' '},num2str(dims(Idim)),{'D '},methodStr(Imet));
%         
%             %Sort by test type
%             xVec = (xRange).^(1/dims(Idim));
%             yVec = time_avg(testMask);
%             yVec_max = time_max(testMask);
%             yVec_min = time_min(testMask);
%             %Sort by proc Label
%             xVec = xVec(procMask);
%             yVec = yVec(procMask);
%             yVec_max = yVec_max(procMask);
%             yVec_min = yVec_min(procMask);
% 
%             %Ordering
%             [xVec,I]=sort(xVec);
%             yVec = yVec(I);
%             yVec_max = yVec_max(I);
%             yVec_min = yVec_min(I);
%             npVec = nb_procs(I);
% 
%             %Lower Bound For Plotting
%             lBoundPlot = 1;
%             hBoundPlot = numel(yVec);
%             %if(dims(Idim) == 3 && strcmp(methodStr(Imet), 'FFT-l'))
%             if(strcmp(Legend(end), 'WEAK 3D FFT-l'))
%                 if (trial_processors)
%                     Legend{end} = ['P=', num2str(p_label(i))];
%                 end
%                 %lBoundPlot = 2;
%             end
%             
%             if(strcmp(Legend(end), 'COMP 3D FFT-g'))
%                 if (trial_processors)
%                     Legend{end} = ['P=', num2str(p_label(i))];
%                 end
%                 %lBoundPlot = 2;
%             end
%             
%             if(strcmp(Legend(end), 'WEAK 2D FFT-g'))
%                 if (trial_processors)
%                     Legend{end} = ['P=', num2str(p_label(i))];
%                 end
%                 %hBoundPlot = 8;
%             end
% 
%             xVec = xVec(lBoundPlot:hBoundPlot);
%             yVec = yVec(lBoundPlot:hBoundPlot);
%             yVec_min = yVec_min(lBoundPlot:hBoundPlot);
%             yVec_max = yVec_max(lBoundPlot:hBoundPlot);
%             npVec = npVec(lBoundPlot:hBoundPlot);
% 
%             if(numel(xVec)<1)
%                 continue
%             end
% 
%             %Normalization
%             %yVec = yVec/yVec(1);
% 
%             pointTags = cell(numel(yVec),1);
%             for i = 1:numel(yVec)
%                 if(strcmp(methodStr{Imet}(end), 'l'))
%                     %pointTags{i} = sprintf('\n\n%d', npVec(i));
%                 else
%                     %pointTags{i} = sprintf('%d\n\n', npVec(i));
%                 end
%             end
% 
%             %Ploting
%             if(plot_MaxMin == 1)
%                 grey = [0.4,0.4,0.4];
%                 polygonX = [xVec; xVec(end:-1:1)];
%                 polygonY = [yVec_max; yVec_min(end:-1:1)];
%                 polygonZ = 0*polygonY -0.1*Imet;
%                 %polygon = fill(polygonX, polygonY, grey);
%                 polygon =  patch(polygonX, polygonY, polygonZ, grey);
%                 set(polygon,'facealpha',.35);
%                 set(get(get(polygon,'Annotation'),'LegendInformation'),...
%                        'IconDisplayStyle','off'); % Exclude line from legend
%             end
% 
%             plot(xVec, yVec, thisLine, 'Color', thisColor,'LineWidth', lWidth,...
%                   'MarkerSize',mSize, 'MarkerEdgeColor','k',...
%                   'MarkerFaceColor',thisFill);
% 
%             baseText = Legend(end);
% 
% 
%             %y_max = plot(xVec, yVec_max, thisLine, 'Color', 'r','LineWidth', lWidth-1);
%             %Legend(end+1) = strcat(baseText, ' MAX');
%             %set(get(get(y_max,'Annotation'),'LegendInformation'),...
%             %       'IconDisplayStyle','off'); % Exclude line from legend
% 
%             %y_min = plot(xVec, yVec_min, thisLine, 'Color', 'b','LineWidth', lWidth-1);
%             %Legend(end+1) = strcat(baseText, ' MIN');
%             %set(get(get(y_min,'Annotation'),'LegendInformation'),...
%             %       'IconDisplayStyle','off'); % Exclude line from legend
% 
% 
% 
% 
%                %markers = scatter(xVec,yVec,mArea, thisMarker,'MarkerEdgeColor',thisMLine,...
%             %      'MarkerFaceColor',thisFill,...
%             %      'LineWidth',1.5);
%             %set(get(get(markers,'Annotation'),'LegendInformation'),...
%             %       'IconDisplayStyle','off'); % Exclude line from legend
%             %plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
%             text(xVec,yVec,pointTags,'HorizontalAlignment','center', 'FontSize', 15);
%         end
%         
%         xlabel('$L/\ell_c$', 'FontSize', 20);
%         %xlabel('Number of processors', 'FontSize', 20);
%         if(time_Type == Wall_Time)
%             ylabel('Wall Time [s]', 'FontSize', 20);
%         elseif (time_Type == CPU_Time)
%             ylabel('CPU Time [s]', 'FontSize', 20);
%         end
%     end
% end
% 
% 
%         
% grid('on')
% box('on')
% set(gca,'xscale','log')
% set(gca,'yscale','log', 'FontSize',15)
% 
% xl = xlim;
% yl = ylim;
% lWidth2 = 2;
% switch testType
%     case 'WEAK'
%         
%     case 'COMP'
%         ylim([0.01 yl(2)]);
%         yl = ylim; 
%         
%         %xRef = xRef.^(1/3);
%         
%         Legend(end+1) = {'O(NLogN)'}; 
%         xRef = linspace(xl(1), xl(2), 100);
%         yRef = xRef(:).*log(xRef(:));        
%         yRef(:) = (yRef(:)/yRef(1));
%         
%         plot(xRef.^(1/3), yRef, 'k:', 'LineWidth', lWidth2);
% 
%         Legend(end+1) = {'O(N)'};
%         yRef = xRef(:);
%         yRef(:) = (yRef(:)/yRef(1));
%         plot(xRef.^(1/3), yRef, 'k--', 'LineWidth', lWidth2);
%         
%         Legend(end+1) = {'O(N^2)'};
%         yRef = (xRef(:)).^2;
%         yRef(:) = (yRef(:)/yRef(1));
%         plot(xRef.^(1/3), yRef, 'k-.', 'LineWidth', lWidth2);
% 
% end
% 
% legend(Legend,'Location','southeast','FontSize',15)
% cd(baseFolder);
% 
% % 
% % %% Ploting data
% % fig = figure(1);
% % hold on
% % hold all
% % 
% % switch testType
% %     case 'C'
% %         legendInfo = cell(2*length(methodChar)*length(dims),1);
% %         legendMask = true(2*length(methodChar)*length(dims),1);
% %     case 'W'
% %         legendInfo = cell(length(methodChar)*length(dims)+1,1);
% %         legendMask = true(length(methodChar)*length(dims)+1,1);
% %         legendInfo{end} = 'Reference';
% %     case 'S'
% %         
% % end 
% % 
% % for i = 1:length(methodChar)
% %     
% %     curMet=methodNb(i);
% %     curIndep = false;
% %     if(strcmp(methodChar{i}(end),'i'))
% %         curIndep = true;
% %     end
% %     
% % 
% %     for j = 1:length(dims)
% %         curDim = dims(j);
% %         
% %         dimText = [num2str(curDim),'D'];
% %         
% %         %Filtering
% %         criteria = (nDim == curDim &...
% %             method == curMet &...
% %             strcmp(testTypeVec, testType)&...
% %             independent == curIndep);        
% %         if(sum(criteria) == 0)
% %             legendMask((i-1)*(2*length(dims))+(2*j-1)) = false;
% %             legendMask((i-1)*(2*length(dims))+(2*j)) = false;
% %             continue
% %         end
% %         
% %         switch testType
% %             case 'C'
% %                 %Legend info
% %                 %legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodChar{i},' ',dimText];
% %                 %legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodChar{i},' ',dimText, ' theoretical'];
% %                 legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodBN{curMet},' ',dimText];
% %                 legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodBN{curMet},' ',dimText, ' theoretical'];
% %                 
% %                 xVec = L_powD(criteria);
% %                 yVec = time(criteria);
% %                 yVec2 = complexity(criteria);
% %                 
% %                 %Ordering
% %                 [xVec,I]=sort(xVec);
% %                 yVec = yVec(I);
% %                 yVec2 = yVec2(I);
% %                 
% %                 %Lower Bound For Plotting
% %                 lBoundPlot = [1, 1, 1];
% %                 if(methodChar{i} == 'S')
% %                     lBoundPlot = [4, 4, 1];
% %                 end
% %                 if(methodChar{i} == 'R')
% %                     lBoundPlot = [4, 4, 1];
% %                 end
% %                 if(methodChar{i} == 'F')
% %                     lBoundPlot = [4, 13, 1];
% %                     %uniVec = 1:2:numel(yVec);
% %                     %yVec2 = yVec2./uniVec;
% %                     
% %                 end
% %                 xVec = xVec(lBoundPlot(curDim):end);
% %                 yVec = yVec(lBoundPlot(curDim):end);
% %                 yVec2 = yVec2(lBoundPlot(curDim):end);
% %                 
% %                 if(methodChar{i} == 'F')
% %                     %yVec2 = yVec;
% %                     
% %                 end
% %                 
% %                 if(length(xVec)<1)
% %                     continue
% %                 end
% %                 
% %                 %Normalization
% %                 %xVec = xVec/xVec(1);
% %                 yVec = yVec/yVec(1);
% %                 yVec2 = yVec2/yVec2(1);
% %                 
% %                 %Ploting
% %                 plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
% %                 plot(xVec, yVec2, '-', 'LineWidth', lWidth);
% %                 xlabel('(L/l_c)^d', 'FontSize', 20);
% %                 ylabel('CPU Time', 'FontSize', 20)
% %                 
% %             case 'W'
% %                 if(curIndep)
% %                     legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText,' with localization'];
% %                 else
% %                     legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText];
% %                 end
% %                 %legendInfo{2*dim}   = [dimText, ' theoretical'];
% %                 
% %                 %Filtering
% %                 xVec = nb_procs(criteria);
% %                 yVec = time(criteria)./(xVec);
% %                 
% %                 %Ordering
% %                 [xVec,I]=sort(xVec);
% %                 yVec = yVec(I);
% %                 toto = nNodes(criteria);
% %                 toto = toto(I);
% %                 %yVec2 = yVec2(I);
% %                 
% %                 %Lower Bound For Plotting
% %                 lBoundPlot = [1, 1, 1];
% %                 if(curMet == SHINOZUKA)
% %                     lBoundPlot = [1, 1, 2];
% %                     %xVec = xVec./2;
% %                     
% %                 elseif(curMet == FFT)
% %                      lBoundPlot = [1, 1, 2];   
% %                 end
% %                 
% %                 xVec = xVec(lBoundPlot(curDim):end);
% %                 yVec = yVec(lBoundPlot(curDim):end);
% %                 %yVec2 = yVec2(lBoundPlot(dim):end);
% %                 
% %                 %if(curMet == FFT && curIndep == true)
% %                 %    yVec(6:end) = yVec(5)*yVec(6:end)/yVec(6);
% %                 %    xVec = xVec(1:end-1);
% %                 %    yVec = yVec(1:end-1);
% %                 %end
% %                 
% %                 if(length(xVec)<1)
% %                     continue
% %                 end
% %                 
% %                 %Normalization
% %                 %xVec = xVec/xVec(1);
% %                 yVec = yVec/yVec(1);
% %                 %yVec2 = yVec2/yVec2(1);
% %                 
% %                 %Ploting
% %                 plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
% %                 %plot(xVec, yVec2, '-');
% %                 xlabel('Number of processors', 'FontSize', 20);
% %                 ylabel('Normalized Wall Time', 'FontSize', 20)
% %             case 'S'
% %         
% %         end
% % 
% %     end
% %  end
% % 
% % switch testType
% %     case 'W'
% %         legendInfo{end} = 'Reference';
% %         plot(xVec, ones(length(xVec),1), '-', 'MarkerSize',10, 'LineWidth', lWidth);
% %         ylim([0.1 1000])
% %         xlim([1 max(xVec)])
% % end
% % 
% % grid('on')
% % box('on')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log', 'FontSize',15)
% % legend(legendInfo(legendMask),'Location','northwest','FontSize',15)
% % % 
% % % saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');
% % 
% % hold off