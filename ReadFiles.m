clc
clear
close all

%% USER
%baseFolder = '/home/lcp/Desktop/RF_Matlab';
%baseFolder = '/mssmat2/home/paludo/Desktop/RF_Matlab';
baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';

%methodChar---------------
%'S' for Shinozuka
%'R' for Randomization
%'I' for Isotropic
%'F', for FFT
%For all the methods put aditional i for independent

%testType-----------------
%'C' for Complexity
%'W' for Weak Scaling
%'S' for Strong Scaling

%methodChar = {'S','R','I'}; 
%methodChar = {'S','Si'};

%methodChar = {'S', 'Si', 'R', 'Ri', 'I','Ii'};
%searchFolder = 'Current/WEAK_Test/NEW';

dims       = [2];
methodChar = {'F','Fi'};

testType = 'W'; %'C' for Complexity, 'W' for Weak Scaling, 'S' for Strong Scaling
searchFolder = 'TESTS_RF_V3/WEAK';

%searchFolder = 'Occigen/AutoTest';

% methodChar = {'S', 'R', 'I'};
% searchFolder = 'WEAK_512';

% methodChar = {'S','R','I'}
% searchFolder = 'Current';
% methodChar = {'I','Ii'};
% searchFolder = 'Current/ISO_Test';
%searchFolder = 'NEWEST_BACKUP';
%searchFolder = 'GOOD_RESULTS_BACKUP';

%---------------------------------

%% CONSTANTS AND INITIALIZATION

%Constants
lWidth = 3; %Line width
nMethods = 4;
nTests = 3;
sizeBTVec = 9;
kAdjust = 1;
periodMult = 1.1;
fileName = 'Sample_Info.h5';
resultsFolder = '/results/';

linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));

%Allocation
methodText = cell(nMethods,1);
methodBN = cell(nMethods,1);
testTypeText = cell(2*nTests,1);
testTypeBN = cell(nTests,1);

%Methods
ISOTROPIC = 1;
methodText{ISOTROPIC} = 'ISO';
methodBN{ISOTROPIC} = 'Isotropic Spectral Method';

SHINOZUKA = 2;
methodText{SHINOZUKA} = 'SHINO'; %'SHINO', 'RANDO' or 'ISO'
methodBN{SHINOZUKA} = 'Spectral Method';

RANDOMIZATION = 3;
methodText{RANDOMIZATION} = 'RANDO';
methodBN{RANDOMIZATION} = 'Randomization';

FFT = 4;
methodText{FFT} = 'FFT';
methodBN{FFT} = 'FFT';

%Tests
COMP = 1;
testTypeBN{COMP} = 'Complexity';
testTypeText{COMP} = '/COMP/';
testTypeText{nTests+COMP} = '/COMP-i/';

STRG = 2;
testTypeBN{STRG} = 'Strong Scaling';
testTypeText{STRG} = '/STRONG/';
testTypeText{nTests+STRG} = '/STRONG-i/';

WEAK = 3;
testTypeBN{WEAK} = 'Weak Scaling';
testTypeText{WEAK} = '/WEAK/';
testTypeText{nTests+WEAK} = '/WEAK-i/';

%% READING FILES

cd([searchFolder])
pathList = subdir(fileName); %Gives all the paths that have a singleGen file
cd(baseFolder);
vecSize = size(pathList,1);

%Allocation
xMin  = cell(vecSize,1);
xMax  = cell(vecSize,1);
nFields = cell(vecSize,1);
procExtent = cell(vecSize,1);
xStep  = cell(vecSize,1);
corrL = cell(vecSize,1);
nFields = cell(vecSize,1);
BT_min  = zeros(vecSize,sizeBTVec);
BT_max  = zeros(vecSize,sizeBTVec);
BT_avg  = zeros(vecSize,sizeBTVec);
BT_stdDev  = zeros(vecSize,sizeBTVec);
nNodes = zeros(vecSize,1);
localizationLevel = zeros(vecSize,1);
xNTotal_Loc = zeros(vecSize,1);
kNTotal_Loc = zeros(vecSize,1);
nb_procs = zeros(vecSize,1);
locLevel = cell(vecSize,1);
nDim = zeros(vecSize,1);
corrMod = zeros(vecSize,1);
method = zeros(vecSize,1);
testMask = false(vecSize,1);
procExtent  = cell(vecSize,1);
testTypeVec = cell(vecSize,1);
nInputsOnFile = zeros(vecSize,1);
independent = false(vecSize,1);

for i = 1 : vecSize
    info = hdf5info(pathList(i).name);
    nAttrib = size(info.GroupHierarchy.Attributes,2);
    nDSet  = size(info.GroupHierarchy.Datasets,2);
    
    %Reading Attributes
    for j = 1 : nAttrib
        if (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'nb_procs'))
            nb_procs(i) = info.GroupHierarchy.Attributes(j).Value;
        elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'nDim'))
            nDim(i) = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Nmc'))
%             Nmc(i) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'method'))
             method(i) = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seedStart'))
%             seedStart(i) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'corrMod'))
             corrMod(i) = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'margiFirst'))
%             margiFirst(i) = info.GroupHierarchy.Attributes(j).Value;
%         %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_CPU_Time'))
%         %    gen_CPU_Time(i) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'gen_WALL_Time'))
             gen_WALL_Time(i) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'nFields'))
             nFields{i} = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'procExtent'))
             procExtent{i} = info.GroupHierarchy.Attributes(j).Value;
%         %elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'seed'))
%         %    seed = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xMaxGlob'))
%             xMaxGlob{i} = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xMinGlob'))
%             xMinGlob{i} = info.GroupHierarchy.Attributes(j).Value;
%         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'xStep'))
%             xStep{i} = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'BT_min'))
             BT_min(i,:) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'BT_max'))
             BT_max(i,:) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'BT_avg'))
             BT_avg(i,:) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'BT_stdDev'))
             BT_stdDev(i,:) = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'corrL'))
             corrL{i} = info.GroupHierarchy.Attributes(j).Value;
         elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'overlap'))
             overlap{i} = info.GroupHierarchy.Attributes(j).Value;
        end
    end
    
    %Finding Test Type
    for test=1:size(testTypeText,1)
        if(~isempty(findstr(pathList(i).name, testTypeText{test})))
            testTypeVec{i} = testTypeText{test}(2);
            break;
        end   
    end
    
end


%% Cleaning vectors from repeated and empty files cases

labels = (1:vecSize)';
time = sum(BT_avg,2);

for file = 1:vecSize
    
    %Treating empty files
    %if(nInputsOnFile(file) == 0)
    %    continue
    %end
    
    %Testing if this test is not redundant
    if(labels(file) ~= file)
        continue
    else
        testMask(file) = true;
    end
    
    n = strfind(pathList(file).name, resultsFolder) - 1;
    
    str1 = pathList(file).name;
    
    for file2 = 1:vecSize
        if(file == file2)
            continue
        end
        
        str2 = pathList(file2).name;
        
        if(strncmpi(str1,str2,n))
            %'They Are In The Same Folder, average them'
            labels(file2) = file;
        end
    end
end

%Taking the average time
var_time = accumarray(labels, time, [], @(x) var(x,1));
time = accumarray(labels, time, [], @(x) mean(x,1));
 
 
% Cleaning vectors from repeated cases
nDim = nDim(testMask);
% xMin  = xMin(testMask);
% xMax  = xMax(testMask);
% xStep = xStep(testMask);
% corrL = corrL(testMask);
time  = time(testMask);
var_time = var_time(testMask);
nb_procs = nb_procs(testMask);
method = method(testMask);
corrMod = corrMod(testMask);
% kNTotal_Loc = kNTotal_Loc(testMask);
% xNTotal_Loc = xNTotal_Loc(testMask);
% usedPathList = pathList(testMask);
% totalSize = sum(testMask);
% testTypeVec = testTypeVec(testMask);
% nInputsOnFile = nInputsOnFile(testMask);
% nNodes = nNodes(testMask);
% var_coef = var_time./time;
% locLevel = locLevel(testMask);
% procExtent  = procExtent(testMask);
nFields = nFields(testMask);
procExtent = procExtent(testMask);
% independent = independent(testMask);
% 
%% Processing data


%Creating methodNb (Mapping from the Char Names to a Number)
methodNb = zeros(length(methodChar),1);
for i = 1:length(methodChar)
    if(strcmp(methodChar{i}(1), 'R'))
        methodNb(i) = RANDOMIZATION;
    elseif(strcmp(methodChar{i}(1), 'S')) 
        methodNb(i) = SHINOZUKA;
    elseif(strcmp(methodChar{i}(1), 'I')) 
        methodNb(i) = ISOTROPIC;
    elseif(strcmp(methodChar{i}(1), 'F')) 
        methodNb(i) = FFT;
    end
end
% 
% %Creating L_powD Vector
% L_powD = zeros(totalSize,1);
% for i = 1:size(xMax,1)
%     L_powD(i) = sqrt(sum(((xMax{i}-xMin{i})./corrL{i}).^2));
%     L_powD(i) = prod(((xMax{i}-xMin{i})./corrL{i}));
% end
% 
% %Creating complexity Vectors
% complexity = zeros(totalSize,1);
% for i = 1:size(xMax,1)
%     complexity(i) = theoretical_complexity(xStep{i}, (xMax{i}-xMin{i})./corrL{i}, corrMod(i), method(i), kAdjust, periodMult);
% end
% 
%Creating independent Vector
independent = true(size(nFields,1),1);
for i = 1:size(nFields,1)
    if(prod(nFields{i}) == 1)
        independent(i) = false;
    end
end


%% Ploting data
fig = figure(1);
hold on
hold all

switch testType
    case 'C'
        legendInfo = cell(2*length(methodChar)*length(dims),1);
        legendMask = true(2*length(methodChar)*length(dims),1);
    case 'W'
        legendInfo = cell(length(methodChar)*length(dims)+1,1);
        legendMask = true(length(methodChar)*length(dims)+1,1);
        legendInfo{end} = 'Reference';
    case 'S'
        
end

for i = 1:length(methodChar)
    
    curMet=methodNb(i);
    curIndep = false;
    if(strcmp(methodChar{i}(end),'i'))
        curIndep = true;
    end
    

    for j = 1:length(dims)
        curDim = dims(j);
        
        dimText = [num2str(curDim),'D'];
        
        %Filtering
        criteria = (nDim == curDim &...
            method == curMet &...
            strcmp(testTypeVec, testType)&...
            independent == curIndep);        
        if(sum(criteria) == 0)
            legendMask((i-1)*(2*length(dims))+(2*j-1)) = false;
            legendMask((i-1)*(2*length(dims))+(2*j)) = false;
            continue
        end
        
        switch testType
            case 'C'
                %Legend info
                %legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodChar{i},' ',dimText];
                %legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodChar{i},' ',dimText, ' theoretical'];
                legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodBN{curMet},' ',dimText];
                legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodBN{curMet},' ',dimText, ' theoretical'];
                
                xVec = L_powD(criteria);
                yVec = time(criteria);
                yVec2 = complexity(criteria);
                
                %Ordering
                [xVec,I]=sort(xVec);
                yVec = yVec(I);
                yVec2 = yVec2(I);
                
                %Lower Bound For Plotting
                lBoundPlot = [1, 1, 1];
                if(methodChar{i} == 'S')
                    lBoundPlot = [4, 4, 1];
                end
                if(methodChar{i} == 'R')
                    lBoundPlot = [4, 4, 1];
                end
                if(methodChar{i} == 'F')
                    lBoundPlot = [4, 13, 1];
                    %uniVec = 1:2:numel(yVec);
                    %yVec2 = yVec2./uniVec;
                    
                end
                xVec = xVec(lBoundPlot(curDim):end);
                yVec = yVec(lBoundPlot(curDim):end);
                yVec2 = yVec2(lBoundPlot(curDim):end);
                
                if(methodChar{i} == 'F')
                    %yVec2 = yVec;
                    
                end
                
                if(length(xVec)<1)
                    continue
                end
                
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
                if(curIndep)
                    legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText,' with localization'];
                else
                    legendInfo{(i-1)*(length(dims))+(j)} = [methodBN{curMet},' ',dimText];
                end
                %legendInfo{2*dim}   = [dimText, ' theoretical'];
                
                %Filtering
                xVec = nb_procs(criteria);
                yVec = time(criteria)./(xVec);
                
                %Ordering
                [xVec,I]=sort(xVec);
                yVec = yVec(I);
                toto = nNodes(criteria);
                toto = toto(I);
                %yVec2 = yVec2(I);
                
                %Lower Bound For Plotting
                lBoundPlot = [1, 1, 1];
                if(curMet == SHINOZUKA)
                    lBoundPlot = [1, 1, 2];
                    %xVec = xVec./2;
                    
                elseif(curMet == FFT)
                     lBoundPlot = [1, 1, 2];   
                end
                
                xVec = xVec(lBoundPlot(curDim):end);
                yVec = yVec(lBoundPlot(curDim):end);
                %yVec2 = yVec2(lBoundPlot(dim):end);
                
                %if(curMet == FFT && curIndep == true)
                %    yVec(6:end) = yVec(5)*yVec(6:end)/yVec(6);
                %    xVec = xVec(1:end-1);
                %    yVec = yVec(1:end-1);
                %end
                
                if(length(xVec)<1)
                    continue
                end
                
                %Normalization
                %xVec = xVec/xVec(1);
                yVec = yVec/yVec(1);
                %yVec2 = yVec2/yVec2(1);
                
                %Ploting
                plot(xVec, yVec, '--^', 'MarkerSize',10, 'LineWidth', lWidth);
                %plot(xVec, yVec2, '-');
                xlabel('Number of processors', 'FontSize', 20);
                ylabel('Normalized Wall Time', 'FontSize', 20)
            case 'S'
        
        end

    end
 end

switch testType
    case 'W'
        legendInfo{end} = 'Reference';
        plot(xVec, ones(length(xVec),1), '-', 'MarkerSize',10, 'LineWidth', lWidth);
        ylim([0.1 1000])
        xlim([1 max(xVec)])
end

grid('on')
box('on')
set(gca,'xscale','log')
set(gca,'yscale','log', 'FontSize',15)
legend(legendInfo(legendMask),'Location','northwest','FontSize',15)
% 
% saveas(fig,[testTypeBN, '_', methodBN,'_L'],'epsc');

hold off