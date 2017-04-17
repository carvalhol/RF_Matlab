clc
close all
clear all

set(0,'defaulttextinterpreter','latex')

partition_type = 'T'; %T = tri, C = cos

alpha = 1;
l_c = 1;
l_weigth = 2;
sz = 100;
x = alpha*l_c*(0:sz)/sz;

Orange = [1 0.3 0];
Blue   = [0 0 1];

figure(1);
switch(partition_type)
    case 'T'
        phi = 1 - x/alpha/l_c; phi(x/alpha/l_c>1) = 0; %Triangular
        color=Orange;
        hold all
        plot(x, phi,'color',color, 'LineWidth',l_weigth);
        plot(x, phi(end:-1:1),'--','color',color, 'LineWidth',l_weigth);
        text(0.32,0.8,'$~~\psi(x^*) = \displaystyle{1 - x^*}$', 'FontSize', 33);
        text(0.1,0.5,'TRI', 'FontSize', 33,'Color',color, 'FontWeight', 'bold');
        fig_name='Phi_tri';
    case 'C'
        phi = 1/2+(1/2*(cos(x*pi/alpha/l_c))); phi(x/alpha/l_c>1) = 0; %Cosinus
        color=Blue;
        hold all
        plot(x, phi,'color',color, 'LineWidth',l_weigth);
        plot(x, phi(end:-1:1),'--','color',color, 'LineWidth',l_weigth);
        text(0.32,0.8,'$~~\psi(x^*) = \displaystyle{\frac{1}{2} + \frac{\cos(\pi x^*)}{2}}$', 'FontSize', 33);
        text(0.1,0.5,'COS', 'FontSize', 33,'Color',color, 'FontWeight', 'bold');
        fig_name='Phi_cos';
end

%xlabel('$x$', 'FontSize', 25);
%ylabel('Wall Time [s]', 'FontSize', 20);
set(gca,'xtick',[0, 1],'ytick',[0, 1], 'LineWidth', l_weigth, 'FontSize', 33);
set(1, 'Position', [0, 0, 650, 250]);
rule_fig(1)
saveas(1,fig_name,'epsc');