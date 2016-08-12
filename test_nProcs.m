close all
clear all
clc

set(0,'defaulttextinterpreter','latex')

size = 11;
L_loc = 32;
L_2 = (2*ones(1,size-log2(L_loc))).^(log2(L_loc):size-1);
L = L_loc:.01:2^(size-1);
lc =1;
d = 3;

N = (L/lc).^d;
N_2 = (L_2/lc).^d;
N_loc = (L_loc/lc).^d;
%Nprocs = ceil((N.*log(N))/(N_loc.*log(N_loc)));
Nprocs = ((N.*log(N))/(N_loc.*log(N_loc)));
Nprocs_loc = ((N_2.*log(N_2))/(N_loc.*log(N_loc)));

fontSize = 25;
markerSize = 15;
lineSize = 5;

figure(1)
hold on
plot(L, Nprocs,'LineWidth', lineSize);
plot(L, (L/L(1)).^d,'LineWidth', lineSize);
plot(L_2, Nprocs_loc,'+','MarkerSize', markerSize, 'MarkerFaceColor','b',...
    'LineWidth', lineSize, 'MarkerEdgeColor','k');
grid on
box on
xlabel('$L/l_c$', 'FontSize', fontSize);
ylabel('Number of processors', 'FontSize', fontSize)
set(gca, 'xtick', (2.^(log2(L_loc):(size-1)))) % set ticks at 1,2,4,8,...
set(gca,'xscale','log', 'FontSize',fontSize-5)
set(gca,'yscale','log', 'FontSize',fontSize-5)
xlim([L_loc 2^size]); 

widthFig=500;
heightFig=400;

set(1, 'Position', [100, 100, widthFig, heightFig]);
rule_fig(1)
saveas(1,'test_nProcs','epsc');