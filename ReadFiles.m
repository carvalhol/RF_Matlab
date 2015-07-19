clc
clear
close all

%% USER
baseFolder = '/home/lcp/Desktop/RF_Matlab';
% baseFolder = '/mssmat2/home/paludo/Desktop/RF_Matlab';

testType = 'W'; %'C' for Complexity, 'W' for Weak Scaling, 'S' for Strong Scaling
%methodChar = {'S','R','I'}; %'S' for Shinozuka, 'R' for Randomization, 'I' for Isotropic, put aditional i for independent
%methodChar = {'S','Si'};
dims       = [3];

% methodChar = {'S', 'Si', 'R', 'Ri', 'I','Ii'};
% searchFolder = 'Current/WEAK_Test/NEW';

methodChar = {'S', 'Si', 'R', 'Ri'};
searchFolder = 'little_WEAK';

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
nMethods = 3;
nTests = 3;
sizeTimeVec = 3;
kAdjust = 1;
periodMult = 1.1;
fileName = 'singleGen';
resultsFolder = 'results';

linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));

%Allocation
methodText = cell(nMethods,1);
methodBN = cell(nMethods,1);
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

%Tests
COMP = 1;
testTypeBN{COMP} = 'Complexity';

STRG = 2;
testTypeBN{STRG} = 'Strong Scaling';

WEAK = 3;
testTypeBN{WEAK} = 'Weak Scaling';

%% READING FILES

cd([searchFolder])
pathList = subdir(fileName); %Gives all the paths that have a singleGen file
cd(baseFolder);
vecSize = size(pathList,1);

%Allocation
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
xMin_Loc  = cell(vecSize,1);
xMax_Loc  = cell(vecSize,1);
indep  = false(vecSize,1);
testTypeVec = cell(vecSize,1);

for file = 1:vecSize
    %pathList(file).name
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
        elseif(strcmp(line, ' --independent-----------------------'))
            data = fgetl(fid);
            indep(file) = strcmp(data, ' T');                        
        end        
    end
    %Finding Test Type
    n = findstr(pathList(file).name, resultsFolder) - 2;
    str1 = pathList(file).name;
    testTypeVec{file} = str1(n);
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
    
    n = findstr(pathList(file).name, resultsFolder) - 1;
    
    str1 = pathList(file).name;
    
    for file2 = 1:vecSize
        if(file == file2)
            continue
        end
        
        str2 = pathList(file2).name;
        
        if(strncmpi(str1,str2,n))
            %'They Are In The Same Folder, take just one'
            labels(file2) = file;
        end
    end
end

%Taking the average time
time = accumarray(labels, time, [], @(x) mean(x,1));

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
kNTotal_Loc = kNTotal_Loc(testMask);
xNTotal_Loc = xNTotal_Loc(testMask);
usedPathList = pathList(testMask);
totalSize = sum(testMask);
indep = indep(testMask);
testTypeVec = testTypeVec(testMask);

%% Processing data

%Creating methodNb Vector
methodNbVec = zeros(totalSize,1);
methodNbVec(strcmp(method, ' RANDOMIZATION')) = RANDOMIZATION;
methodNbVec(strcmp(method, ' SHINOZUKA')) = SHINOZUKA;
methodNbVec(strcmp(method, ' ISOTROPIC')) = ISOTROPIC;

%Creating methodNb
methodNb = zeros(length(methodChar),1);
for i = 1:length(methodChar)
    if(strcmp(methodChar{i}(1), 'R'))
        methodNb(i) = RANDOMIZATION;
    elseif(strcmp(methodChar{i}(1), 'S')) 
        methodNb(i) = SHINOZUKA;
    elseif(strcmp(methodChar{i}(1), 'I')) 
        methodNb(i) = ISOTROPIC;
    end
end

%Creating L_powD Vector
L_powD = zeros(totalSize,1);
for i = 1:size(xMax,1)
    L_powD(i) = sqrt(sum(((xMax{i}-xMin{i})./corrL{i}).^2));
    L_powD(i) = prod(((xMax{i}-xMin{i})./corrL{i}));
end

%Creating complexity Vectors
complexity = zeros(totalSize,1);
for i = 1:size(xMax,1)
    complexity(i) = theoretical_complexity(xStep{i}, (xMax{i}-xMin{i})./corrL{i}, corrMod{i}, kAdjust, periodMult);
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
            methodNbVec == curMet &...
            strcmp(testTypeVec, testType) &...
            indep == curIndep);        
        if(sum(criteria) == 0)
            legendMask((i-1)*(2*length(dims))+(2*j-1)) = false;
            legendMask((i-1)*(2*length(dims))+(2*j)) = false;
            continue
        end
        
        switch testType
            case 'C'
                %Legend info
                legendInfo{(i-1)*(2*length(dims))+(2*j-1)} = [methodChar{i},' ',dimText];
                legendInfo{(i-1)*(2*length(dims))+(2*j)}   = [methodChar{i},' ',dimText, ' theoretical'];
                
                
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
                    lBoundPlot = [8, 6, 2];
                end
                if(methodChar{i} == 'R')
                    lBoundPlot = [1, 1, 1];
                end
                xVec = xVec(lBoundPlot(curDim):end);
                yVec = yVec(lBoundPlot(curDim):end);
                yVec2 = yVec2(lBoundPlot(curDim):end);
                
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
                %yVec2 = yVec2(I);
                
                %Lower Bound For Plotting
                lBoundPlot = [1, 1, 1];
                xVec = xVec(lBoundPlot(curDim):end);
                yVec = yVec(lBoundPlot(curDim):end);
                %yVec2 = yVec2(lBoundPlot(dim):end);
                
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