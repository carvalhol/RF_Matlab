clc
close all
clear all

sz = 3;
y_lim = [-1.2 1.2];
l_weigth = 3;
l_weigth_2 = 2;
FontSize_y = 15;
FontSize_x = 18;
TitleSize = 25;
widthFig = 270;
heightFig = sz *100;

kPoints = linspace(0.5*pi,2*pi,sz);
xPoints = linspace(0,10,100);
randField = zeros(numel(xPoints, 1));
phi_hat = 2*pi*rand(1, sz);
Sk = exp(-kPoints.^2/(4*pi));

%plot(xPoints, cos(kPoints(1)*xPoints));
%hold on
%hold all

set(0,'defaulttextinterpreter','latex')

cc=lines(12);

%g(\xv) = \sum_{\nv \leq \Nv} \sqrt{2\mathcal{S}(\kv_{\nv})|\Delta\kv|}\cos{( \kv_{\nv}\cdot \xv + \hat{\phi}_{\nv})}

figure(1); 
hold on;
count = 0;
for n=1:sz
    subplot(sz,1,n+count);
    plot(xPoints, cos(kPoints(n)*xPoints),'color',cc(n,:), 'LineWidth',l_weigth);
    if(n == 1)
        title('$\cos(\mathbf{k}_\mathbf{n}\cdot\mathbf{x})$','FontSize',TitleSize)
    end
    ylim(y_lim)
    xlabel('$\mathbf{x}$','FontSize',FontSize_x)
    ylabel({['$\mathbf{k_',num2str(n),'}=',num2str(kPoints(n)/pi),' \pi$'],' ','Amplitude'},'FontSize',FontSize_y)
    box on;
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
end

set(1, 'Position', [100, 100, widthFig, heightFig]);
rule_fig(1)
saveas(1,'RF_Step_1-Cos(kx)','epsc');



figure(2); 
hold on;
count = 0;
for n=1:sz
    subplot(sz,1,n+count);
    plot(xPoints, cos(kPoints(n)*xPoints+phi_hat(n)),'color',cc(n,:), 'LineWidth',l_weigth);
    if(n == 1)
        title('$\cos(\mathbf{k}_\mathbf{n}\cdot\mathbf{x}+\hat{\phi}_\mathbf{n})$','FontSize',TitleSize)
    end
    ylim(y_lim)
    xlabel('$\mathbf{x}$','FontSize',FontSize_x)
    %ylabel({['$k=',num2str(kPoints(n)/pi),' \pi$'],' ','Amplitude'},'FontSize',FontSize_y)
    box on;
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
end

set(2, 'Position', [100, 100, widthFig, heightFig]);
rule_fig(2)
saveas(2,'RF_Step_2-Cos(kx+phi)','epsc');

figure(3); 
hold on;
count = 0;
for n=1:sz
    subplot(sz,1,n+count);
    plot(xPoints, sqrt(Sk(n))*cos(kPoints(n)*xPoints+phi_hat(n)),'color',cc(n,:), 'LineWidth',l_weigth);
    if(n == 1)
        title('$\hat{\mathbf{R}}^{1/2}_\mathbf{n}\cos(\mathbf{k}_\mathbf{n}\cdot\mathbf{x}+\hat{\phi}_\mathbf{n})$','FontSize',TitleSize)
    end
    ylim(y_lim)
    xlabel('$\mathbf{x}$','FontSize',FontSize_x)
    %ylabel({['$k=',num2str(kPoints(n)/pi),' \pi$'],' ','Amplitude'},'FontSize',FontSize_y)
    box on;
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
    randField=randField+sqrt(Sk(n))*cos(kPoints(n)*xPoints+phi_hat(n));
end

set(3, 'Position', [100, 100, widthFig, heightFig]);
rule_fig(3)
saveas(3,'RF_Step_3-Sk_Cos(kx+phi)','epsc');

figure(4);
plot(xPoints, randField,'color','r', 'LineWidth',l_weigth);
title(' ')
%ylim(y_lim)
xlabel('$\mathbf{x}$','FontSize',FontSize_x)
ylabel({['$\mu(\mathbf{x})$']},'FontSize',FontSize_y)
box on;
set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
randField=randField+sqrt(Sk(n))*cos(kPoints(n)*xPoints+phi_hat(n));

set(4, 'Position', [100, 100, widthFig, heightFig/sz]);
rule_fig(4)
saveas(4,'RF_Step_4-randField','epsc');

close all


