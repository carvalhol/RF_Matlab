%
% _monitor parameter_
%
xs = [];
ys = [];
zs = [];

%GRID Premiere couche (Fonctions de coh?rence)
x = -200:10:200;
y = -50:10:50;
z = 0;

[xv, yv, zv] = meshgrid(x,y,z);
xs = [xs, reshape(xv,[1, numel(xv)])];
ys = [ys, reshape(yv,[1, numel(yv)])];
zs = [zs, reshape(zv,[1, numel(zv)])];

[ys, I] = sort(ys);
xs = xs(I);
zs = zs(I);

%PLAN XZ (Fonctions de transfert)
x = -200:10:200;
y = 0;
z = [-300;-150];

[xv, yv, zv] = meshgrid(x,y,z);
xs = [xs, reshape(xv,[1, numel(xv)])];
ys = [ys, reshape(yv,[1, numel(yv)])];
zs = [zs, reshape(zv,[1, numel(zv)])];

%PLAN YZ (Fonctions de transfert)
x = 0;
y = -200:10:200;
z = [-300;-150];

[xv, yv, zv] = meshgrid(x,y,z);
xs = [xs, reshape(xv,[1, numel(xv)])];
ys = [ys, reshape(yv,[1, numel(yv)])];
zs = [zs, reshape(zv,[1, numel(zv)])];

%SOURCE
xv = 0;
yv = 0;
zv = -1400;
xs = [xs, reshape(xv,[1, numel(xv)])];
ys = [ys, reshape(yv,[1, numel(yv)])];
zs = [zs, reshape(zv,[1, numel(zv)])];


monitor.name 			= 'SC3D';
monitor.type 			= 'points';
monitor.fname 			= 'stations.txt';
monitor.period 			= 41;
monitor.x               = xs;
monitor.y               = ys;
monitor.z               = zs;
monitor.nm              = numel(monitor.z);


SC3D_generate_stationstxt(monitor);