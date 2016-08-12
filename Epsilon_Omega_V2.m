clc
clear
close all

%% USER
set(0,'defaulttextinterpreter','latex')
baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab';
fileName = 'Sample_Info.h5';

searchFolder = '/Users/carvalhol/Desktop/ERROR_TEST/2D/FFT-g';
imageName = '/Users/carvalhol/Desktop/GITs/Bouvier2014_GIT/figures/FFT-g';
nFolders = 3;

% searchFolder = '/Users/carvalhol/Desktop/ERROR_TEST/2D/RAN-g';
% imageName = '/Users/carvalhol/Desktop/GITs/Bouvier2014_GIT/figures/RAN-g';
% nFolders = 2;

% searchFolder = '/Users/carvalhol/Desktop/ERROR_TEST/2D/SHI-g';
% imageName = '/Users/carvalhol/Desktop/GITs/Bouvier2014_GIT/figures/SHI-g';
% nFolders = 2;

makerStyle = {'--'; '-'; ':'};
mean_1 = zeros(nFolders,1);
mean_2 = zeros(nFolders,1);
legendInfo_1 = cell(nFolders,1);
legendInfo_2 = cell(nFolders,1);
widthFig = 500;
heightFig = 400;

%for i = 3:3
for i = 1:nFolders
    searchFolderIt = strcat(searchFolder,'/', sprintf('%03d',i));
    cd([searchFolderIt])
    
    %Saving (for better performance)
    values_name = 'pathList';
    fullPath = [pwd,'/',values_name,'.mat'];
    if(exist(fullPath, 'file') == 2)
        load(values_name)
    else
        pathList = subdir(fileName); %Gives all the paths that have a Sample_Info.h5 file
        save(values_name,values_name)      
    end
    
    vecSize = size(pathList,1);

    %Saving (for better performance)
    values_name = 'Sk';
    fullPath = [pwd,'/',values_name,'.mat'];
    if(exist(fullPath, 'file') == 2)
        
        load(values_name)
        
    else
        
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
        end        
        save(values_name,'Sk_1', 'Sk_2')      
    end 
    
    
    for file = 1:vecSize
        
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
%     edges = [-5 -4 -2 -1 -0.5 0 0.5 1 2 4 5];
%     N = histcounts(Distance_L2_Sk1,edges);
    [h_1, centers_1] = hist(Distance_L2_Sk1,nLevels); %Utiliser histc
    [h_2, centers_2] = hist(Distance_L2_Sk2,nLevels);
    dx_1 = centers_1(2) - centers_1(1);
    dx_2 = centers_2(2) - centers_2(1);
    %hist_1 = histogram(Distance_L2_Sk1,nLevels);
    %h_1 = hist_1.Values;
    %centers_1 = hist_1.BinEdges(1:end-1)+hist_1.BinWidth/2;
    %area_1=sum(h_1)*hist_1.BinWidth;
    width_1=centers_1(end)-centers_1(1);
    area_1= dx_1*vecSize;
    %area_1= vecSize;
    h_1=h_1/area_1;
    mean_1(i)=mean(Distance_L2_Sk1);
    
    %hist_2 = histogram(Distance_L2_Sk2,nLevels);
    %h_2 = hist_2.Values;
    %centers_2 = hist_2.BinEdges(1:end-1)+hist_2.BinWidth/2;
    %area_2=sum(h_2)*hist_2.BinWidth;
    width_2=centers_2(end)-centers_2(1);
    area_2= dx_2*vecSize;
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
            %yl = ylim;
            %yl = [0 2000000];
            %line([mean_1(j) mean_1(j)],yl,'LineStyle', lineStyle,'LineWidth', lWidth*j, 'Color','b')
            legendInfo_1{j} = ['$N_k 10^', num2str(j-1),'$'];
            hold off
            
            figure(2)
            hold on
            %yl = ylim;
            %yl = [0 2000000];
            %line([mean_2(j) mean_2(j)],yl,'LineStyle', lineStyle, 'LineWidth', lWidth*j, 'Color','b')
            legendInfo_2{j} = ['$N_k 10^', num2str(j-1),'$'];
            hold off
        end 
        
        cd(baseFolder);
        
        figure(1)
        hold on
        set(gca,'FontSize',20)
        xlabel('L2 Distance', 'FontSize', 25);
        ylabel('Frequency', 'FontSize', 25);       
        ylim([0,20]);
        xlim([0,1.5]);
        box on;
        grid on;
        leg=legend(legendInfo_1,'Location','northeast','FontSize',25,'box','on');
        set(leg,'Interpreter','latex');
        set(1, 'Position', [100, 100, widthFig, heightFig]);
        rule_fig(1)
        saveas(1,[imageName,'_1'],'epsc');
        hold off

        figure(2)
        hold on
        set(gca,'FontSize',20)
        xlabel('L2 Distance', 'FontSize', 25);
        ylabel('Frequency', 'FontSize', 25);
        ylim([0,20]);
        xlim([0,1.5]);
        box on;
        grid on;
        leg=legend(legendInfo_2,'Location','northeast','FontSize',25,'box','on');
        set(leg,'Interpreter','latex');
        set(2, 'Position', [100, 100, widthFig, heightFig]);
        rule_fig(2)
        saveas(2,[imageName,'_2'],'epsc');
        hold off
    end
    
end