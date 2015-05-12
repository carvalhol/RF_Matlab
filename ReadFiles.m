clc
clear
close all

%% USER
testType = 'C'; %'C' for Complexity, 'W' for Weak Scaling, 'S' for Strong Scaling
%nDim = 2;
methodChar = 'S'; %'S' for Shinozuka, 'R' for Randomization, 'I' for Isotropic
sep = '/'; % '/' for linux and Mac, '\' for Windows
equalEdg = false; %Only plot domains where all the Dimensions are of the same size
deformDomainSize = true;
%imposeSearchFolder = false;
imposeSearchFolder = true;
%searchFolder = 'WEAK'; %This argument will be ignored if "imposeSearchFolder = false"
%searchFolder = 'COMP';
searchFolder = 'GOOD_RESULTS_BACKUP/COMP';
lWidth = 3;
%---------------------------------

% CONSTANTS AND INITIALIZATION
RANDOMIZATION = 3;
SHINOZUKA = 2;
ISOTROPIC = 1;

%mask = cell(3,1);
%domainSize = cell(3,1);
%totalTerms = cell(3,1);
%genTime = cell(3,1);
%writeTime = cell(3,1);
sizeTimeVec = 3;
kAdjust = 1;
periodMult = 1.1;

fileName = 'singleGen';
resultsFolder = 'results';

switch testType
    case 'C'
        %rootFolder = 'COMP';
        testTypeBN = 'Complexity';
    case 'S'
        %rootFolder = 'STRONG';
        testTypeBN = 'Strong Scaling';
    case 'W'
        %rootFolder = 'WEAK';
        testTypeBN = 'Weak Scaling';
end

if(~imposeSearchFolder)
    switch testType
        case 'C'
            searchFolder = 'COMP';
        case 'S'
            searchFolder = 'STRONG';
        case 'W'
            searchFolder = 'WEAK';
    end
end

switch methodChar
    case 'S'
        methodText = 'SHINO'; %'SHINO', 'RANDO' or 'ISO'
        methodBN = 'Spectral Method';
        methodNb = SHINOZUKA;
    case 'R'
        methodText = 'RANDO';
        methodBN = 'Randomization';
        methodNb = RANDOMIZATION;
    case 'I'
        methodText = 'ISO';
        methodBN = 'Isotropic Spectral Method';
        methodNb = ISOTROPIC;
end



%% READING FILES


%cd(['GOOD_RESULTS_BACKUP/',searchFolder])
cd([searchFolder])
pathList = subdir(fileName); %Gives all the paths that have a singleGen file
cd('/mssmat2/home/paludo/Desktop/RF_Matlab');
vecSize = size(pathList,1);

    %     dim_xMin  = zeros(nDim, nIter);
    %     dim_xMax  = zeros(nDim, nIter);
    %     dim_corrL = zeros(nDim, nIter);
    %     dim_time  = zeros(3, nIter);

xMin  = cell(vecSize,1);
xMax  = cell(vecSize,1);
xStep  = cell(vecSize,1);
corrL = cell(vecSize,1);
timeVec  = zeros(vecSize,sizeTimeVec);
sum_xNTotal = zeros(vecSize,1);
sum_kNTotal = zeros(vecSize,1);
xNTotal_Loc = zeros(vecSize,1);
kNTotal_Loc = zeros(vecSize,1);
nb_procs = zeros(vecSize,1);
nDim = zeros(vecSize,1);
method = cell(vecSize,1);
testMask = false(vecSize,1);
complexity = zeros(vecSize,1);
xMin_Loc  = cell(vecSize,1);
xMax_Loc  = cell(vecSize,1);

for file = 1:vecSize
    fid = fopen(pathList(file).name);
    line = fgetl(fid);
    while(line ~= -1)
        line = fgetl(fid);
        if(strcmp(line, ' --nDim-----------------------'))
            data = fgetl(fid);
            nDim(file) = str2num(data);
        elseif(strcmp(line, ' --xMinGlob-----------------------'))
            data = fgetl(fid);
            xMin{file} = str2num(data);
        elseif(strcmp(line, ' --xMaxGlob-----------------------'))
            data = fgetl(fid);
            xMax{file} = str2num(data);
        elseif(strcmp(line, ' --xStep-----------------------'))
            data = fgetl(fid);
            xStep{file} = str2num(data);
        elseif(strcmp(line, ' --xMin_Loc-----------------------'))
            data = fgetl(fid);
            xMin_Loc{file} = str2num(data);
        elseif(strcmp(line, ' --xMax_Loc-----------------------'))
            data = fgetl(fid);
            xMax_Loc{file} = str2num(data);
        elseif(strcmp(line, ' --corrL-----------------------'))
            data = fgetl(fid);
            corrL{file} = str2num(data);
        elseif(strcmp(line, ' --timeVec-----------------------'))            
            for i = 1 : sizeTimeVec
                data = fgetl(fid);
                timeVec(file,i) = str2num(data);
            end
        elseif(strcmp(line, ' --sum_xNTotal-----------------------'))
            data = fgetl(fid);
            sum_xNTotal(file) = str2num(data);
        elseif(strcmp(line, ' --sum_kNTotal-----------------------'))
            data = fgetl(fid);
            sum_kNTotal(file) = str2num(data);
        elseif(strcmp(line, ' --xNTotal_Loc-----------------------'))
            data = fgetl(fid);
            xNTotal_Loc(file) = str2num(data);
        elseif(strcmp(line, ' --kNTotal_Loc-----------------------'))
            data = fgetl(fid);
            kNTotal_Loc(file) = str2num(data);
        elseif(strcmp(line, ' --nb_procs-----------------------'))
            data = fgetl(fid);
            nb_procs(file) = str2num(data);             
        elseif(strcmp(line, ' --method-----------------------'))
            data = fgetl(fid);
            method{file} = data; 
        elseif(strcmp(line, ' --corrMod-----------------------'))
            data = fgetl(fid);
            corrMod{file} = data;    
        end        
    end
end

%% Cleaning vectors from repeated cases

labels = (1:vecSize)';
time = timeVec(:,2) - timeVec(:,1);

for file = 1:vecSize
    %Testing if this test is not redundant
    if(labels(file) ~= file)
        continue
    else
        testMask(file) = true;
    end
%     %Testing if it is the good kind of Test
%     if (size(findstr(pathList(file).name, rootFolder),1) > 0) 
%         testMask(file) = true;
%     else
%         continue
%     end
    
    n = findstr(pathList(file).name, resultsFolder) - 1;
    for file2 = 1:vecSize
        if(file == file2)
            continue
        end
        str1 = pathList(file).name;
        str2 = pathList(file2).name;
        
        if(strncmpi(str1,str2,n))
            %'They Are In The Same Folder, take just one'
            labels(file2) = file;
        end
    end
end

%Taking the average time
time = accumarray(labels,time,[],@(x) mean(x,1));

% Cleaning vectors from repeated cases
nDim = nDim(testMask);
xMin  = xMin(testMask);
xMax  = xMax(testMask);
xStep = xStep(testMask);
xMin_Loc = xMin_Loc(testMask);
xMax_Loc = xMax_Loc(testMask);
corrL = corrL(testMask);
time  = time(testMask);
sum_xNTotal = sum_xNTotal(testMask);
sum_kNTotal = sum_kNTotal(testMask);
nb_procs = nb_procs(testMask);
method = method(testMask);
corrMod = corrMod(testMask);
complexity = complexity(testMask);
kNTotal_Loc = kNTotal_Loc(testMask);
xNTotal_Loc = xNTotal_Loc(testMask);
usedPathList = pathList(testMask);

%% Processing data
indicator1 = (fix(sum_xNTotal./xNTotal_Loc) == nb_procs);
indicator2 = (sum_kNTotal./kNTotal_Loc == nb_procs);

%Creating L_powD Vector
L_powD = zeros(size(nDim,1),1);

for i = 1:size(xMax,1)
    %L_powD(i) = sqrt(sum(((xMax{i}-xMin{i})./corrL{i}).^2));
    L_powD(i) = prod(((xMax{i}-xMin{i})./corrL{i}));
end

%Creating complexity Vectors
nTerms = sum_xNTotal.* (sum_kNTotal./nb_procs);
for i = 1:size(xMax,1)
    complexity(i) = theoretical_complexity(xStep{i}, (xMax{i}-xMin{i})./corrL{i}, corrMod{i}, kAdjust, periodMult);
end
%Creating methodNb Vector
methodNbVec = zeros(size(method,1),1);
for i = 1:size(methodNbVec,1)
    if(strcmp(method{i}, ' RANDOMIZATION'))
        methodNbVec(i) = RANDOMIZATION;
    elseif(strcmp(method{i}, ' SHINOZUKA'))
        methodNbVec(i) = SHINOZUKA;
    elseif(strcmp(method{i}, ' ISOTROPIC'))
        methodNbVec(i) = ISOTROPIC;
    end
end

%% Ploting data

fig = figure(1);
hold on
hold all

for dim = 1:3
    
    dimText = [num2str(dim),'D'];    
       
    switch testType
        case 'C'
            %Legend info
            if(dim == 1)
                legendInfo = cell(6,1);
            end
            legendInfo{2*dim-1} = [dimText];
            legendInfo{2*dim}   = [dimText, ' theoretical'];
    
            %Filtering
            xVec = L_powD(nDim == dim & methodNbVec == methodNb);            
            yVec = time(nDim == dim & methodNbVec == methodNb);
            yVec2 = complexity(nDim == dim & methodNbVec == methodNb);
            
            %Ordering
            [xVec,I]=sort(xVec);
            yVec = yVec(I);
            yVec2 = yVec2(I);
            
            %Lower Bound For Plotting
            lBoundPlot = [1, 1, 1];
            if(methodChar == 'S')
                lBoundPlot = [8, 5, 2];
            end
            if(methodChar == 'R')
                lBoundPlot = [7, 5, 1];
            end
            xVec = xVec(lBoundPlot(dim):end);
            yVec = yVec(lBoundPlot(dim):end);
            yVec2 = yVec2(lBoundPlot(dim):end);
            
            %Normalization
            %xVec = xVec/xVec(1); 
            yVec = yVec/yVec(1);
            yVec2 = yVec2/yVec2(1);
            
            %Ploting
            plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
            plot(xVec, yVec2, '-', 'LineWidth', lWidth);
            xlabel('(L/l_c)^d', 'FontSize', 20);
            ylabel('CPU Time', 'FontSize', 20) 
            
        case 'W'
            %Legend info
            if(dim == 1)
                legendInfo = cell(4,1);
            end
            legendInfo{dim} = [dimText];
            %legendInfo{2*dim}   = [dimText, ' theoretical'];
    
            %Filtering
            xVec = nb_procs(nDim == dim & methodNbVec == methodNb);            
            yVec = time(nDim == dim & methodNbVec == methodNb)./xVec;
            %yVec2 = complexity(nDim == dim & methodNbVec == methodNb)./xVec;
            
            %Ordering
            [xVec,I]=sort(xVec);
            yVec = yVec(I);
            %yVec2 = yVec2(I);
            
            %Lower Bound For Plotting
            lBoundPlot = [2, 2, 2];
            xVec = xVec(lBoundPlot(dim):end);
            yVec = yVec(lBoundPlot(dim):end);
            %yVec2 = yVec2(lBoundPlot(dim):end);
            
            %Normalization
            %xVec = xVec/xVec(1); 
            yVec = yVec/yVec(1);
            %yVec2 = yVec2/yVec2(1);
            
            %Ploting
            plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
            %plot(xVec, yVec2, '-');
            xlabel('Number of processors', 'FontSize', 20);
            ylabel('CPU Time / Number of processors', 'FontSize', 20) 
        case 'S'
                       
    end
    
    %existMethods = unique(method);
    
    
    
end

if(testType == 'W')
    limH = xLim;
    limV = yLim;
    xLim([0, limH(2)]);
    yLim([0.1, limV(2)]);
    %axis([0 inf 0.1 1000])
    limH = xLim;
    plot([1, limH(2)], [1,1], 'LineWidth', lWidth);
    legendInfo{4} = ['Reference'];
end

grid('on')
box('on')
set(gca,'xscale','log')
set(gca,'yscale','log', 'FontSize',15)
legend(legendInfo,'Location','southeast','FontSize',15)

saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');

hold off

%% View a specific item

%View
item = 1;

if((item <= size(nDim,1)) && item > 0)
    '--------------------------------------'
    view_item = item
    view_pathList = pathList(item).name
    view_nDim = nDim(item)
    view_xMin = xMin{item}
    view_xMax = xMax{item}
    view_xStep = xStep{item}
    view_xMin_Loc = xMin_Loc{item}
    view_xMax_Loc = xMax_Loc{item}
    view_corrL = corrL{item}
    view_time = time(item)
    view_sum_xNTotal = sum_xNTotal(item)
    view_sum_kNTotal = sum_kNTotal(item)
    view_nb_procs = nb_procs(item)
    view_method = method(item)
    view_corrMod = corrMod(item)
    view_complexity = complexity(item)
    view_kNTotal_Loc = kNTotal_Loc(item)
    view_xNTotal_Loc = xNTotal_Loc(item)
else
    'Item out of range'
end

%% TRASH

% if(deformDomainSize)
%     saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');
% else
%     saveas(fig,[testTypeBN, '_', methodBN,'_LpowD'],'epsc');
% end

%for a = 1:3
    
    %     switch testType
    %         case 'C'
    %             rootFolder = 'COMP';
    %             testTypeBN = 'Complexity';
    %         case 'S'
    %             rootFolder = 'STRONG';
    %             testTypeBN = 'Strong Scaling';
    %         case 'W'
    %             rootFolder = 'WEAK';
    %             testTypeBN = 'Weak Scaling';
    %     end
    %
    %     switch method
    %         case 'S'
    %             methodText = 'SHINO'; %'SHINO', 'RANDO' or 'ISO'
    %             methodBN = 'Spectral Method';
    %         case 'R'
    %             methodText = 'RANDO';
    %             methodBN = 'Randomization';
    %         case 'I'
    %             methodText = 'ISO';
    %             methodBN = 'Isotropic Spectral Method';
    %     end
    %
    %     %Defining Names
    %     dim = [num2str(nDim),'D']; %Ex: "3D"
    %     dimFolder = [dim,'_',testType]; %Ex: "3D_C"
    %     testFolder = [methodText,'_',dim,'_',testType]; %Ex: "SHINO_3D_C"
    %     path = [rootFolder,sep,dimFolder,sep,testFolder,sep]; %Ex: "COMP/3D_C/SHINO_3D_C"
    %     resFolder = 'results';
    %
    %
    
    %
    %     index = dir (path);
    %     nIter = 0;
    %
    %     %Finding Max Iteration Index
    %     for i = 1:size(index,1)
    %         %sprintf(index(i).name(1:3))
    %         %sprintf(['it',num2str(it)])
    %         if(size(index(i).name(:),1) < 4)
    %         elseif(strcmp(index(i).name(1:2), 'it'))
    %             for undIn = 4:10
    %                 if(strcmp(index(i).name(undIn), '_'))
    %                     iter = str2num(index(i).name(3:undIn-1)); %#ok<*ST2NM>
    %                     if(iter > nIter)
    %                         nIter = iter;
    %                     end
    %                     break
    %                 end
    %             end
    %         end
    %     end
    %
    %     %Initializing
    %     dim_xMin  = zeros(nDim, nIter);
    %     dim_xMax  = zeros(nDim, nIter);
    %     dim_corrL = zeros(nDim, nIter);
    %     dim_time  = zeros(3, nIter);
    %     dim_xNTotal = zeros(1, nIter);
    %     dim_kNTotal = zeros(1, nIter);
    %     dim_nb_procs = zeros(1, nIter);
    %     mask = true(1, nIter);
    %     domainSize = zeros(1, nIter);
    %     totalTerms = zeros(1, nIter);
    %     genTime = zeros(1, nIter);
    %     writeTime = zeros(1, nIter);
    %
    %     [testTypeBN,' ',methodBN,' ',dim]
    %
    %     for i = 1:size(index,1)
    %
    %         %Finding file
    %         if(size(index(i).name,2) < 3)
    %             %'Iteration not found'
    %             continue
    %         else
    %             %'inside else'
    %             for it = 1:size(index,1)
    %                 %sprintf(index(i).name(1:3))
    %                 %sprintf(['it',num2str(it)])
    %                 if(strcmp(index(i).name(1:2+size(num2str(it),2)), ['it',num2str(it)]))
    %                     %'Iteration found in pos '
    %                     itFolder = index(i).name;
    %                     for undIn = 4:10
    %                         if(strcmp(index(i).name(undIn), '_'))
    %                             iter = str2num(index(i).name(3:undIn-1)); %#ok<*ST2NM>
    %                             break
    %                         end
    %                     end
    %                     break
    %                 end
    %             end
    %
    %             itPath = [rootFolder,sep,dimFolder,sep,testFolder,sep,itFolder,sep];
    %
    %             %Check for results folder existence
    %             if(exist([itPath,resFolder], 'dir') ~= 7)
    %                 mask(iter) = false;
    %                 continue
    %             end
    %
    %             tempIndex = dir ([itPath,resFolder]);
    %
    %             filePathFound = false;
    %             for it = 1:size(tempIndex,1)
    %                 if(size(tempIndex(it).name,2) < 3)
    %                     %'This is not the good folder'
    %                     continue
    %                 else
    %                     szN = size(tempIndex(it).name,2);
    %                     if(strcmp(tempIndex(it).name(szN-3:szN), '_res'))
    %                         %sprintf(tempIndex(it).name)
    %                         filePath = [itPath,resFolder, sep,tempIndex(it).name,sep,'singleGen'];
    %                         filePathFound = true;
    %                     end
    %                 end
    %             end
    %
    %             if(~filePathFound)
    %                 mask(iter) = false;
    %             end
    %         end
    %
    %
    %
    %         if(~mask(iter))
    %             continue
    %         end
    %
    %         if(exist(filePath, 'file') ~= 2)
    %             mask(iter) = false;
    %             continue
    %         end
    %
    %         %Reading File
    %         fid = fopen(filePath);
    %         line = fgetl(fid);
    %         while(line ~= -1)
    %             line = fgetl(fid);
    %             if(strcmp(line, ' --xMinGlob-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_xMin(:,iter) = str2num(data);
    %                 %             for j = 1 : nDim
    %                 %                 data = fgetl(fid);
    %                 %                 all_xMin(j,iter) = str2num(data);
    %                 %             end
    %             end
    %             if(strcmp(line, ' --xMaxGlob-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_xMax(:,iter) = str2num(data);
    %                 %             for j = 1 : nDim
    %                 %                 data = fgetl(fid);
    %                 %                 all_xMax(j,iter) = str2num(data);
    %                 %             end
    %             end
    %             if(strcmp(line, ' --corrL-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_corrL(:,iter) = str2num(data);
    %                 %             for j = 1 : nDim
    %                 %                 data = fgetl(fid);
    %                 %                 all_corrL(j,iter) = str2num(data);
    %                 %             end
    %             end
    %             if(strcmp(line, ' --timeVec-----------------------'))
    %                 %             data = fgetl(fid);
    %                 %             all_time(:,iter) = str2num(data);
    %                 for j = 1 : size(dim_time,1)
    %                     data = fgetl(fid);
    %                     dim_time(j,iter) = str2num(data);
    %                 end
    %             end
    %
    %             if(strcmp(line, ' --sum_xNTotal-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_xNTotal(1,iter) = str2num(data);
    %             end
    %
    %             if(strcmp(line, ' --sum_kNTotal-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_kNTotal(1,iter) = str2num(data);
    %             end
    %
    %             if(strcmp(line, ' --nb_procs-----------------------'))
    %                 data = fgetl(fid);
    %                 dim_nb_procs(1,iter) = str2num(data);
    %             end
    %         end
    %     end
    %
    %     dim_xMin;
    %     dim_xMax;
    %     dim_nb_procs;
    %     dim_kNTotal;
    %     dim_xNTotal;
    %
    %     if(equalEdg)
    %         mask(mask) = ((prod(dim_xMax(:,mask) - dim_xMin(:,mask),1) ...
    %             == (dim_xMax(1,mask)- dim_xMin(1,mask)).^nDim)) & mask(mask);
    %     end
    %
    %     domainSize(mask) = prod(dim_xMax(:,mask) - dim_xMin(:,mask),1);
    %     totalTerms(mask) = dim_kNTotal(mask).*dim_xNTotal(mask)./dim_nb_procs(mask); %TEST ./dim_nb_procs(mask)
    %     genTime(mask) = dim_time(2,mask) - dim_time(1,mask);
    %     writeTime(mask) = dim_time(3,(mask)) - dim_time(2,(mask));
    %
    %     domainNorm = domainSize(mask);
    %     domainNorm = domainNorm(1);
    %
    %     %% Stocking Data
    %     all_xMin{nDim} = dim_xMin;
    %     all_xMax{nDim}  = dim_xMax;
    %     all_corrL{nDim} = dim_corrL;
    %     all_time{nDim}  = dim_time;
    %     all_xNTotal{nDim} = dim_xNTotal;
    %     all_kNTotal{nDim} = dim_kNTotal;
    %     all_nb_procs{nDim} = dim_nb_procs;
    %     all_domainSize{nDim} = domainSize;
    %     all_totalTerms{nDim} = totalTerms;
    %     all_genTime{nDim} = genTime;
    %     all_writeTime{nDim} = writeTime;
    %     all_mask{nDim} = mask;
    %
%end

% %% PLOTING FIGURES
% fig = figure(2);
% hold on
% hold all
%
% for nDim = 1:3
%
%     dim = [num2str(nDim),'D'];
%     legendInfo{2*nDim-1} = [dim];
%     legendInfo{2*nDim}   = [dim, ' theoretical'];
%     if(deformDomainSize)
%         defDim = nDim;
%     else
%         defDim = 1;
%     end
%
%     if(testType == 'C')
%         if(nDim == 1)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%         if(nDim == 2)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%         if(nDim == 3)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%     elseif(testType == 'W')
%         if(nDim == 1)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%         if(nDim == 2)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%         if(nDim == 3)
%             maskInf = 0;
%             all_mask{nDim}(1:maskInf) = false;
%
%         end
%     end
%
%     timeNorm = all_genTime{nDim}(all_mask{nDim});
%     timeNorm = timeNorm(1);
%     termsNorm = all_totalTerms{nDim}(all_mask{nDim});
%     termsNorm = termsNorm(1);
%
%     domainVec = (all_domainSize{nDim}(all_mask{nDim})/domainNorm).^(1/defDim);
%     timeVec = all_genTime{nDim}(all_mask{nDim})./...
%               (timeNorm*all_nb_procs{nDim}(all_mask{nDim})); %Time per proc
%     termsVec = all_totalTerms{nDim}(all_mask{nDim})./...
%              (termsNorm*all_nb_procs{nDim}(all_mask{nDim})); %Terms per proc
%
%     if(testType == 'C')
%
%     elseif(testType == 'W')
%         timeVec = timeVec./(domainVec.^(defDim));
%         termsVec = termsVec./(domainVec.^(defDim));
%
%         domainVec = all_nb_procs{nDim}(all_mask{nDim});
%     else
%
%     end
%
%     plot(domainVec, timeVec, '--^', 'MarkerSize',10);
%     plot(domainVec, termsVec, '-');
% end
%
% if(deformDomainSize)
%     xlabel('L/l_c', 'FontSize', 20);
% else
%     xlabel('(L/l_c)^d', 'FontSize', 20);
% end
%
% ylabel('Normalized Generation Time', 'FontSize', 20)
%
% if(testType == 'W')
%     xlabel('Number of processors', 'FontSize', 20);
%     ylabel('Wall Time / nProc', 'FontSize', 20);
%     ylim([0.5 2]);
% end
%
% grid('on')
% box('on')
%
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% legend(legendInfo,'Location','southeast')
% %legend('boxoff')
% %axis equal
% hold off
%
% if(deformDomainSize)
%     saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');
% else
%     saveas(fig,[testTypeBN, '_', methodBN,'_LpowD'],'epsc');
% end

