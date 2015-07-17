clc
clear
close all

%% FFT 1D

% lc = 1;
% t = 0:0.01:10;
% nIter = 100
% corrMod = 'gaussian'
% fig = figure(1);
% hold on
% hold all
% 
% for i = 1 : nIter
%     y=randomgenerationfieldsFFTMeth(t, 'gaussian', lc);
%     plot(y)
% end
% 
% grid('on')
% box('on')
% hold off

%% FFT 2D

% x = 0:0.1:10;
% y = 0:0.1:10;
% corrMod = 'gaussian';
% corrModx = corrMod;
% corrMody = corrMod;
% delta = 0.01;
% Lx = 1;
% Ly = 1;
% 
% nIter = 1;
% 
% fig = figure(1);
% hold on
% hold all
% 
% for i = 1 : nIter
%     u=ranfomfieldsgeneration2DFFT(x,y,corrModx,corrMody , Lx, Ly, delta);
%     surf(u)
%     shading flat
%     view(2)
%     colorbar
% end
% 
% grid('on')
% box('on')
% hold off

%% FFT 3D

Lmax = 4;

x = 0:0.1:Lmax;
y = 0:0.1:Lmax;
z = 0:0.1:Lmax;
corrMod = 'gaussian';
corrModx = corrMod;
corrMody = corrMod;
corrModz = corrMod;
delta = 0.01;
Lx = 1;
Ly = 1;
Lz = 1;

nIter = 1;

for i = 1 : nIter
    %u=ranfomfieldsgeneration3DFFT(x,y,z,corrModx,corrMody,corrModz,Lx,Ly,Lz,delta);    
    u=Test_3DFFT(x,y,z,corrModx,corrMody,corrModz,Lx,Ly,Lz,delta);        
end

%% Results plot

xMinPlot = 1;
yMinPlot = 1;
zMinPlot = 1;
% xLimPlot = 20;
% yLimPlot = 20;
% zLimPlot = 20;
xLimPlot = 10000;
yLimPlot = 10000;
zLimPlot = 10000;
xStepPlot = 3;
yStepPlot = 3;
zStepPlot = 3;

if(xLimPlot > numel(x))
    xLimPlot = numel(x);
end
if(yLimPlot > numel(y))
    yLimPlot = numel(y);
end
if(zLimPlot > numel(z))
    zLimPlot = numel(z);
end
if(xMinPlot < 1)
    xMinPlot = 1;
end
if(yMinPlot < 1)
    yMinPlot = 1;
end
if(zMinPlot < 1)
    zMinPlot = 1;
end
%uV = u(xMin:xLim, yMin:yLim, zMin:zLim);
% [xV,yV,zV] = meshgrid(x(xMin:xLim),y(yMin:yLim),z(zMin:zLim));

[xV,yV,zV] = meshgrid(...
             x(xMinPlot:xStepPlot:xLimPlot),...
             y(yMinPlot:yStepPlot:yLimPlot),...
             z(zMinPlot:zStepPlot:zLimPlot));

uV = u(...
    xMinPlot:xStepPlot:xLimPlot,...
    yMinPlot:yStepPlot:yLimPlot,...
    zMinPlot:zStepPlot:zLimPlot);

% xV = xV();
% yV = yV();
% zV = zV;

A = figure(1);
%set(A, 'render', 'zbuffer')
scatter3(xV(:), yV(:), zV(:), 5,uV(:))

std(u(:))
mean(u(:))