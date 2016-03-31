clc
clear
close all

%% USER
%baseFolder = '/home/lcp/Desktop/RF_Matlab';
%baseFolder = '/mssmat2/home/paludo/Desktop/RF_Matlab';
baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';

%% User Entries

%methodStr---------------
%'SHI' for Shinozuka
%'RAN' for Randomization
%'ISO' for Isotropic
%'FFT', for FFT
%For all the methods choose -g for global and -l for localization method

% searchFolder = 'NEW_Tests/WEAK_2';
% testType   = 'WEAK'; %'COMP' for Complexity, 'WEAK' for Weak Scaling, 'STRONG' for Strong Scaling
% dims       = [2,3]; %Which dimensions take into account.
% methodStr = {'FFT-g','FFT-l'};

searchFolder = 'NEW_Tests';
testType   = 'COMP'; %'COMP' for Complexity, 'WEAK' for Weak Scaling, 'STRONG' for Strong Scaling
dims       = [2,3]; %Which dimensions take into account.
methodStr = {'FFT-g','FFT-l'};


%---------------------------------

%% CONSTANTS AND INITIALIZATION

%Constants
nMethods = 4;
kAdjust = 1;
fileName = 'Sample_Info.h5';
resultsFolder = '/results/';
sizeBTVec = 9;
nTests = numel(dims)*numel(methodStr);
Legend=cell(0,1);

linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));

%% READING FILES AND PLOTTING DATA

fig = figure(1);
hold on
hold all

for Idim = 1:numel(dims)
    for Imet = 1:numel(methodStr)
        
        cd(baseFolder);
        cd(searchFolder);
        
        cropedPath = strcat('./', ...
                     testType,'/', ...
                     num2str(dims(Idim)),'D/', ...
                     methodStr(Imet));
                 
        Legend(end+1) = strcat(testType,{' '},num2str(dims(Idim)),{'D '},methodStr(Imet));
                 
        cd(cropedPath{1});
        folderList = subdir(fileName);
        vecSize = numel(folderList);

        %Allocation---------------------
        nb_procs = zeros(vecSize,1);
        nDim = zeros(vecSize,1);
        method = zeros(vecSize,1);
        corrMod = zeros(vecSize,1);
        gen_WALL_Time = zeros(vecSize,1);
        nFields  = zeros(vecSize,dims(Idim));
        procExtent = zeros(vecSize,dims(Idim));
        xMaxGlob = zeros(vecSize,dims(Idim));
        xMinGlob = zeros(vecSize,dims(Idim));
        xStep    = zeros(vecSize,dims(Idim));
        BT_min  = zeros(vecSize,sizeBTVec);
        BT_max  = zeros(vecSize,sizeBTVec);
        BT_avg  = zeros(vecSize,sizeBTVec);
        BT_stdDev  = zeros(vecSize,sizeBTVec);
        corrL   = zeros(vecSize,dims(Idim));
        overlap = zeros(vecSize,dims(Idim));
        localizationLevel = zeros(vecSize,1);

        for Ipath = 1:numel(folderList)
            
            folderList(Ipath).name
            info = hdf5info(folderList(Ipath).name);
            nAttrib = size(info.GroupHierarchy.Attributes,2);
            nDSet  = size(info.GroupHierarchy.Datasets,2);
            %file_id = H5F.open(filename);
            %H5F.close(file_id);
            
            %% Reading Attributes
            for Iattr = 1 : nAttrib
                if (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nb_procs'))
                    nb_procs(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nDim'))
                    nDim(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
                    %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Nmc'))
                    %             Nmc(i) = info.GroupHierarchy.Attributes(j).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'method'))
                    method(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
                    %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seedStart'))
                    %             seedStart(i) = info.GroupHierarchy.Attributes(j).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'corrMod'))
                    corrMod(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
                    %         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'margiFirst'))
                    %             margiFirst(i) = info.GroupHierarchy.Attributes(j).Value;
                    %         %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_CPU_Time'))
                    %         %    gen_CPU_Time(i) = info.GroupHierarchy.Attributes(j).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'gen_WALL_Time'))
                    gen_WALL_Time(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'nFields'))
                    nFields(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'procExtent'))
                    procExtent(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                    %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seed'))
                    %    seed = info.GroupHierarchy.Attributes(j).Value;
                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xMaxGlob'))
                     xMaxGlob(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                 elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xMinGlob'))
                     xMinGlob(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'xStep'))
                    xStep(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_min'))
                    BT_min(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_max'))
                    BT_max(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_avg'))
                    BT_avg(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'BT_stdDev'))
                    BT_stdDev(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'corrL'))
                    corrL(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'overlap'))
                    overlap(Ipath,:) = info.GroupHierarchy.Attributes(Iattr).Value;
                elseif (strcmp(info.GroupHierarchy.Attributes(Iattr).Shortname, 'localizationLevel'))
                    localizationLevel(Ipath) = info.GroupHierarchy.Attributes(Iattr).Value;    
                end
            end
        end
        
        %% Cleaning vectors from repeated and empty files cases
        
        testMask = false(vecSize, 1);
        
        labels = (1:vecSize)';
        time = sum(BT_avg(:,1:3),2).*prod(nFields.^repmat(localizationLevel,1,dims(Idim)),2);
        %time = gen_WALL_Time;
        
        for file = 1:vecSize

            %Testing if this test is not redundant
            if(labels(file) ~= file)
                continue
            else
                testMask(file) = true;
            end
            
            n = strfind(folderList(file).name, resultsFolder) - 1;
            
            str1 = folderList(file).name;
            
            for file2 = 1:vecSize
                if(file == file2)
                    continue
                end
                
                str2 = folderList(file2).name;
                
                if(strncmpi(str1,str2,n))
                    %'They Are In The Same Folder, average them'
                    labels(file2) = file;
                end
            end
        end
        
        %Taking the time average and variance
        var_time = accumarray(labels, time, [], @(x) var(x,1));
        time = accumarray(labels, time, [], @(x) mean(x,1));

        %% Plotting Results
        
        lWidth = 5;
        
        %switch testType
            %case 'WEAK'               
                %Filtering
                nPoints = prod(1+int64((xMaxGlob(testMask,:) - xMinGlob(testMask,:))./xStep(testMask,:)),2);
                xRange = prod((xMaxGlob(testMask,:) - xMinGlob(testMask,:))./corrL(testMask,:),2);
                nb_procs = nb_procs(testMask);
                
                xVec = xRange;
                yVec = time(testMask);
                
                
                %Ordering
                [xVec,I]=sort(xVec);
                yVec = yVec(I);
                
                %Lower Bound For Plotting
                lBoundPlot = 1;
                hBoundPlot = numel(yVec);
                %if(dims(Idim) == 3 && strcmp(methodStr(Imet), 'FFT-l'))
                if(strcmp(Legend(end), 'WEAK 3D FFT-l'))
                    lBoundPlot = 2;
                end
                if(strcmp(Legend(end), 'WEAK 2D FFT-g'))
                    hBoundPlot = 8;
                end
                
                xVec = xVec(lBoundPlot:hBoundPlot);
                yVec = yVec(lBoundPlot:hBoundPlot);
                nb_procs = nb_procs(lBoundPlot:hBoundPlot);
                
                if(numel(xVec)<1)
                    continue
                end
                
                %Normalization
                %yVec = yVec/yVec(1);
                
                pointTags = cell(numel(yVec),1);
                for i = 1:numel(yVec)
                    pointTags{i} = sprintf('\n\n%d', nb_procs(i));
                end
                
                %Ploting
                plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
                text(xVec,yVec,pointTags,'HorizontalAlignment','center', 'FontSize', 15);
                xlabel('(L/l_c)^d', 'FontSize', 20);
                %xlabel('Number of processors', 'FontSize', 20);
                ylabel('Wall Time [s]', 'FontSize', 20);
            %case 'COMP'
                
        %end
    end
end

grid('on')
box('on')
set(gca,'xscale','log')
set(gca,'yscale','log', 'FontSize',15)
legend(Legend,'Location','southeast','FontSize',15)
cd(baseFolder);

% 
% %% Ploting data
% fig = figure(1);
% hold on
% hold all
% 
% switch testType
%     case 'C'
%         legendInfo = cell(2*length(methodChar)*length(dims),1);
%         legendMask = true(2*length(methodChar)*length(dims),1);
%     case 'W'
%         legendInfo = cell(length(methodChar)*length(dims)+1,1);
%         legendMask = true(length(methodChar)*length(dims)+1,1);
%         legendInfo{end} = 'Reference';
%     case 'S'
%         
% end 
% 
% for i = 1:length(methodChar)
%     
%     curMet=methodNb(i);
%     curIndep = false;
%     if(strcmp(methodChar{i}(end),'i'))
%         curIndep = true;
%     end
%     
% 
%     for j = 1:length(dims)
%         curDim = dims(j);
%         
%         dimText = [num2str(curDim),'D'];
%         
%         %Filtering
%         criteria = (nDim == curDim &...
%             method == curMet &...
%             strcmp(testTypeVec, testType)&...
%             independent == curIndep);        
%         if(sum(criteria) == 0)
%             legendMask((i-1)*(2*length(dims))+(2*j-1)) = false;
%             legendMask((i-1)*(2*length(dims))+(2*j)) = false;
%             continue
%         end
%         
%         switch testType
%             case 'C'
%                 %Legend info
%                 %legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodChar{i},' ',dimText];
%                 %legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodChar{i},' ',dimText, ' theoretical'];
%                 legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodBN{curMet},' ',dimText];
%                 legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodBN{curMet},' ',dimText, ' theoretical'];
%                 
%                 xVec = L_powD(criteria);
%                 yVec = time(criteria);
%                 yVec2 = complexity(criteria);
%                 
%                 %Ordering
%                 [xVec,I]=sort(xVec);
%                 yVec = yVec(I);
%                 yVec2 = yVec2(I);
%                 
%                 %Lower Bound For Plotting
%                 lBoundPlot = [1, 1, 1];
%                 if(methodChar{i} == 'S')
%                     lBoundPlot = [4, 4, 1];
%                 end
%                 if(methodChar{i} == 'R')
%                     lBoundPlot = [4, 4, 1];
%                 end
%                 if(methodChar{i} == 'F')
%                     lBoundPlot = [4, 13, 1];
%                     %uniVec = 1:2:numel(yVec);
%                     %yVec2 = yVec2./uniVec;
%                     
%                 end
%                 xVec = xVec(lBoundPlot(curDim):end);
%                 yVec = yVec(lBoundPlot(curDim):end);
%                 yVec2 = yVec2(lBoundPlot(curDim):end);
%                 
%                 if(methodChar{i} == 'F')
%                     %yVec2 = yVec;
%                     
%                 end
%                 
%                 if(length(xVec)<1)
%                     continue
%                 end
%                 
%                 %Normalization
%                 %xVec = xVec/xVec(1);
%                 yVec = yVec/yVec(1);
%                 yVec2 = yVec2/yVec2(1);
%                 
%                 %Ploting
%                 plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
%                 plot(xVec, yVec2, '-', 'LineWidth', lWidth);
%                 xlabel('(L/l_c)^d', 'FontSize', 20);
%                 ylabel('CPU Time', 'FontSize', 20)
%                 
%             case 'W'
%                 if(curIndep)
%                     legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText,' with localization'];
%                 else
%                     legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText];
%                 end
%                 %legendInfo{2*dim}   = [dimText, ' theoretical'];
%                 
%                 %Filtering
%                 xVec = nb_procs(criteria);
%                 yVec = time(criteria)./(xVec);
%                 
%                 %Ordering
%                 [xVec,I]=sort(xVec);
%                 yVec = yVec(I);
%                 toto = nNodes(criteria);
%                 toto = toto(I);
%                 %yVec2 = yVec2(I);
%                 
%                 %Lower Bound For Plotting
%                 lBoundPlot = [1, 1, 1];
%                 if(curMet == SHINOZUKA)
%                     lBoundPlot = [1, 1, 2];
%                     %xVec = xVec./2;
%                     
%                 elseif(curMet == FFT)
%                      lBoundPlot = [1, 1, 2];   
%                 end
%                 
%                 xVec = xVec(lBoundPlot(curDim):end);
%                 yVec = yVec(lBoundPlot(curDim):end);
%                 %yVec2 = yVec2(lBoundPlot(dim):end);
%                 
%                 %if(curMet == FFT && curIndep == true)
%                 %    yVec(6:end) = yVec(5)*yVec(6:end)/yVec(6);
%                 %    xVec = xVec(1:end-1);
%                 %    yVec = yVec(1:end-1);
%                 %end
%                 
%                 if(length(xVec)<1)
%                     continue
%                 end
%                 
%                 %Normalization
%                 %xVec = xVec/xVec(1);
%                 yVec = yVec/yVec(1);
%                 %yVec2 = yVec2/yVec2(1);
%                 
%                 %Ploting
%                 plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
%                 %plot(xVec, yVec2, '-');
%                 xlabel('Number of processors', 'FontSize', 20);
%                 ylabel('Normalized Wall Time', 'FontSize', 20)
%             case 'S'
%         
%         end
% 
%     end
%  end
% 
% switch testType
%     case 'W'
%         legendInfo{end} = 'Reference';
%         plot(xVec, ones(length(xVec),1), '-', 'MarkerSize',10, 'LineWidth', lWidth);
%         ylim([0.1 1000])
%         xlim([1 max(xVec)])
% end
% 
% grid('on')
% box('on')
% set(gca,'xscale','log')
% set(gca,'yscale','log', 'FontSize',15)
% legend(legendInfo(legendMask),'Location','northwest','FontSize',15)
% % 
% % saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');
% 
% hold off