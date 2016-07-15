clc
clear
close all

%% USER
%baseFolder = '/home/lcp/Desktop/RF_Matlab';
%baseFolder = '/mssmat2/home/paludo/Desktop/RF_Matlab';
baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';

searchFolder = '/Users/carvalhol/Desktop/ERROR_TEST/2D/SHI-g';
fileName = 'Sample_Info.h5';

nFolders = 2;

makerStyle = {'--'; '-'; ':'};
mean_1 = zeros(nFolders,1);
mean_2 = zeros(nFolders,1);
legendInfo_1 = cell(nFolders,1);
legendInfo_2 = cell(nFolders,1);

%for i = 3:3
for i = 1:nFolders
    searchFolderIt = strcat(searchFolder,'/', sprintf('%03d',i));
    cd([searchFolderIt])
    pathList = subdir(fileName); %Gives all the paths that have a Sample_Info.h5 file
    cd(baseFolder);
    vecSize = size(pathList,1);
    
    Sk_1  = cell(vecSize,1);
    Sk_2  = cell(vecSize,1);
    
    for file = 1:vecSize
        info = hdf5info(pathList(file).name);
        nAttrib = size(info.GroupHierarchy.Attributes,2);
        nDSet  = size(info.GroupHierarchy.Datasets,2);
        
        for j = 1 : nAttrib
            if (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Sk_1'))
                Sk_1{file} = info.GroupHierarchy.Attributes(j).Value;
            elseif (strcmp(info.GroupHierarchy.Attributes(j).Shortname, 'Sk_2'))
                Sk_2{file} = info.GroupHierarchy.Attributes(j).Value;
            end
        end
        
        if(file == 1)
            R_ref1 = zeros(size(Sk_1{1},1),1);
            R_ref2 = zeros(size(Sk_2{1},1),1);
        end
        
        R_ref1 = R_ref1 + Sk_1{file};
        R_ref2 = R_ref2 + Sk_2{file};
    end
    
    pathList(file).name
    R_ref1 = R_ref1/vecSize;
    R_ref2 = R_ref2/vecSize;
    
    Distance_L2_Sk1 = zeros(vecSize, 1);
    Distance_L2_Sk2 = zeros(vecSize, 1);
    
    for sample = 1:vecSize
        Distance_L2_Sk1(sample) = sqrt(sum((Sk_1{sample}-R_ref1).^2)); %Waiting for deltaK multiplication
        Distance_L2_Sk2(sample) = sqrt(sum((Sk_2{sample}-R_ref2).^2));
    end
    
    
    
    nLevels = 20;
    lWidth = 1;
    [h_1, centers_1] = hist(Distance_L2_Sk1,nLevels);
    [h_2, centers_2] = hist(Distance_L2_Sk2,nLevels);
    
    %hist_1 = histogram(Distance_L2_Sk1,nLevels);
    %h_1 = hist_1.Values;
    %centers_1 = hist_1.BinEdges(1:end-1)+hist_1.BinWidth/2;
    %area_1=sum(h_1)*hist_1.BinWidth;
    width_1=centers_1(2)-centers_1(1);
    %area_1= width_1*sum(h_1);
    area_1= vecSize;
    h_1=100*h_1/area_1;
    mean_1(i)=mean(Distance_L2_Sk1);
    
    %hist_2 = histogram(Distance_L2_Sk2,nLevels);
    %h_2 = hist_2.Values;
    %centers_2 = hist_2.BinEdges(1:end-1)+hist_2.BinWidth/2;
    %area_2=sum(h_2)*hist_2.BinWidth;
    width_2=centers_2(2)-centers_2(1);
    area_2= vecSize;
    h_2=h_2/area_2;
    mean_2(i)=mean(Distance_L2_Sk2);
    
    lineStyle = makerStyle{mod(i-1, size(makerStyle,1))+1};
    
    figure(1)
    hold on
    plot(centers_1, h_1, lineStyle, 'MarkerSize',2, 'LineWidth', lWidth*i, 'Color','k');    
    hold off
    
    figure(2)
    hold on
    plot(centers_2, h_2, lineStyle, 'MarkerSize',2, 'LineWidth', lWidth*i, 'Color','k');    
    hold off
    
    %Ploting Average Lines
    if(i == nFolders)
        for j=1:nFolders
            lineStyle = makerStyle{mod(j-1, size(makerStyle,1))+1};
            
            figure(1)
            hold on
            yl = ylim;
            line([mean_1(j) mean_1(j)],yl,'LineStyle', lineStyle,'LineWidth', lWidth*j, 'Color','b')
            legendInfo_1{j} = strcat('Nk_x * 10^', num2str(j-1));
            hold off
            
            figure(2)
            hold on
            yl = ylim;
            line([mean_2(j) mean_2(j)],yl,'LineStyle', lineStyle, 'LineWidth', lWidth*j, 'Color','b')
            legendInfo_2{j} = strcat('Nk_y * 10^', num2str(j-1));
            hold off
        end 
        
        figure(1)
        hold on
        xlabel('L2 Distance', 'FontSize', 20);
        ylabel('Frequency (%)', 'FontSize', 20);
        legend(legendInfo_1,'Location','northeast','FontSize',15)
        hold off

        figure(2)
        hold on
        xlabel('L2 Distance', 'FontSize', 20);
        ylabel('Frequency (%)', 'FontSize', 20);
        legend(legendInfo_2,'Location','northeast','FontSize',15)
        hold off
    end
    
end