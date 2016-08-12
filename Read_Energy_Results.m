close all
clear all
clc

set(0,'defaulttextinterpreter','latex')
Legend=cell(0,1);

baseFolder = '/Users/carvalhol/Desktop/GITs/RF_Matlab/Energy_Results/MOMENT_TEST';
testType  = {'CTE'};
sourceType = {'P'};


figure(1)
hold on

for i = 1:numel(testType)
    cd(baseFolder);
    cd(testType{i});
    list = ls;
    
    for j = 1:numel(sourceType)
        exp = [sourceType{j},'_(\w+)p'];
        [tokens,matches] = regexp(list,'P_(\w+)p','tokens','match');
        int_tokens = zeros(numel(tokens), 1);
        for k = 1:numel(int_tokens)
            int_tokens(k) = str2num(char(tokens{k})); 
        end
        [int_tokens,I]=sort(int_tokens);
        tokens = tokens(I);
        matches = matches(I);
        
        for k = 1:numel(matches)
            cd(matches{k});
            
            fid = fopen('En_P.txt','r');
            formatSpec = [' time=   ' '%f' ' P_En=' '%f'];
            Out_P = fscanf(fid,formatSpec,[2 Inf])';
            fclose(fid);
            Legend(end+1) = strcat('p=',tokens{k},'- P En - Source ',sourceType{j});            
            plot(Out_P(:,1), Out_P(:,2),'LineWidth', 3);                       
            cd('..');
        end
        
        for k = 1:numel(matches)
            cd(matches{k});        
            fid = fopen('En_S.txt','r');
            formatSpec = [' time=   ' '%f' ' S_En=' '%f'];
            Out_S = fscanf(fid,formatSpec,[2 Inf])';
            fclose(fid);
            Legend(end+1) = strcat('p=',tokens{k},'- S En - Source ',sourceType{j});            
            plot(Out_S(:,1), Out_S(:,2),'LineWidth', 3);            
            cd('..');
        end
    end
end

ylabel('Energy', 'FontSize', 20);
xlabel('Time(s)', 'FontSize', 20);
legend(Legend,'Location','southeast','FontSize',15)
