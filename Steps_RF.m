clc
close all
clear all

sz = 4;
y_lim = [-1.2 1.2];
l_weigth = 5;
l_weigth_2 = 5;

kPoints = linspace(pi/sz,pi,sz);
xPoints = linspace(0,10,100);
randField = zeros(numel(xPoints, 1));
phi_hat = 2*pi*rand(1, sz);
Sk = exp(-kPoints.^2/(4*pi));

%plot(xPoints, cos(kPoints(1)*xPoints));
%hold on
%hold all

cc=lines(12);


figure; 
hold on;

count = 0;
for n=1:sz
    subplot(3,sz,n+count)
    plot(xPoints, cos(kPoints(n)*xPoints),'color',cc(n,:), 'LineWidth',l_weigth);
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
    ylim(y_lim)
end

count = sz;
for n=1:sz
    subplot(3,sz,n+count)
    plot(xPoints, cos(kPoints(n)*xPoints+phi_hat(n)),'color',cc(n,:), 'LineWidth',l_weigth);
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
    ylim(y_lim)
end

count = 2*sz;
for n=1:sz
    subplot(3,sz,n+count)
    plot(xPoints, sqrt(Sk(n))*cos(kPoints(n)*xPoints+phi_hat(n)),'color',cc(n,:), 'LineWidth',l_weigth);
    set(gca,'xtick',[],'ytick',[], 'LineWidth', l_weigth_2)
    ylim(y_lim)
    randField = randField + sqrt(Sk(n))*cos(kPoints(n)*xPoints+phi_hat(n));
end

%fig = subplot(3,sz,count:count+sz-1);
%print(fig,'MySavedPlot','-dpng')

%count = 3*sz;
%figure; plot(xPoints, randField);
