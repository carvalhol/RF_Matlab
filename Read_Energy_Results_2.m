close all
clear all
clc

set(0,'defaulttextinterpreter','latex')
Legend=cell(0,1);
Titles=cell(0,1);

baseFolder = '/home/carvalhol/Desktop/Energy';
%baseFolder = '/Users/carvalhol/Desktop/Energy';
testType  = {'CTE'};
sourceType = {'P'};
countFig = 0;



for i = 1:numel(testType)
    cd(baseFolder);
    cd(testType{i});
    list = ls;
    
    for j = 1:numel(sourceType)
        countFig = countFig + 1;
        figure(countFig)
        hold on
        
        Title = strcat('Source: ',sourceType{j}, ', Media: ',testType{i});
        %Titles(end+1) = strcat('Source ',sourceType{j}, ', Media',testType{i});
        title(Title);
        
        exp = [sourceType{j},'_(\w+)p'];
        [tokens,matches] = regexp(list,'P_(\w+)p','tokens','match');
        int_tokens = zeros(numel(tokens), 1);
        for k = 1:numel(int_tokens)
            int_tokens(k) = str2num(char(tokens{k})); 
        end
        [int_tokens,I]=sort(int_tokens);
        tokens = tokens(I);
        matches = matches(I);
        
        for k = numel(matches):-1:1
            cd(matches{k});
            
            fid = fopen('En_P.txt','r');
            formatSpec = [' time=   ' '%f' ' P_En=' '%f'];
            Out_P = fscanf(fid,formatSpec,[2 Inf])';
            fclose(fid);
            Legend(end+1) = strcat('p=',tokens{k},'- P Energy');            
            plot(Out_P(:,1), Out_P(:,2),'LineWidth', 3);                       
        %    cd('..');
        %end
        
        %for k = 1:numel(matches)
        %    cd(matches{k});        
            fid = fopen('En_S.txt','r');
            formatSpec = [' time=   ' '%f' ' S_En=' '%f'];
            Out_S = fscanf(fid,formatSpec,[2 Inf])';
            fclose(fid);
            Legend(end+1) = strcat('p=',tokens{k},'- S Energy');            
            plot(Out_S(:,1), Out_S(:,2),'LineWidth', 3);            
            cd('..');
        end
        
        set(gca, 'yscale','log', 'FontSize',15)
        ylabel('Energy(J)', 'FontSize', 25);
        xlabel('Time(s)', 'FontSize', 25);
        legend(Legend,'Location','northeast','FontSize',20)
        
        hold off
    end
    
end
