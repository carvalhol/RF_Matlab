close all
clear all

L = 50;
x =linspace(0,L,50*L+1);
y = x;
lc = 1;
law = 'gauss';
correl = 'sinc2';
mu = 0.0;
s = 1.0;
Nmc = 1;
zVal = 5;
%modification of seed in space
lowLim=71;
highLim=79;

[xM, yM] = meshgrid(x,y);

k = randomField( law, correl, lc, mu, s, Nmc, x, y, [], [], 3);
k2 = randomField( law, correl, lc, mu, s, Nmc, x, y, [], [], 3, [lowLim; highLim]);
diffField=abs((k2-k));

t = linspace(0,2*pi,1000);
R = 3.2; %Modification spot mark radius
Rcenter = 37*[1 1]; %Modification spot mark radius
cx = Rcenter(1) + R*cos(t);
cy = Rcenter(2) + R*sin(t);
z = 5*cy./cy;

figure('units','normalized','position',[.0 .0 1.0 1.0])
hSurf = surf(xM,yM,k,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
caxis([-zVal,zVal])
hold on; plot3(cx,cy,z,'-k','linewidth',8);
az = 0;
el = 90;
view(az, el);
colormap jet;
%axis off
box on
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
%set(gca,'Visible','off','plotboxaspectratio',[1 1 1]);
set(gca,'plotboxaspectratio',[1 1 1]);
set(gca,'linewidth',8)
print ('-dpng', 'FieldBEFORE')

figure('units','normalized','position',[.0 .0 1.0 1.0])
hSurf = surf(xM,yM,k2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
caxis([-zVal,zVal])
hold on; plot3(cx,cy,z,'-k','linewidth',8);
az = 0;
el = 90;
view(az, el);
colormap jet;
box on
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
set(gca,'plotboxaspectratio',[1 1 1]);
set(gca,'linewidth',8)
print ('-dpng', 'FieldAFTER')


figure('units','normalized','position',[.0 .0 1.0 1.0])
hSurf = surf(xM,yM,diffField,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
hold on; plot3(cx,cy,z,'-k','linewidth',8);
az = 0;
el = 90;
view(az, el);
mpc = colormap('gray');
colormap(mpc(end:-1:1,:));
box on
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
set(gca,'plotboxaspectratio',[1 1 1]);
set(gca,'linewidth',8)
print ('-dpng', 'FieldDIFFERENCE')