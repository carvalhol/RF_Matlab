clc
clear all
close all

law='gauss';
correl='sinc2';
L=10;
mu=0;
nSamples = 1;
s=1;
Nmc=1;
N=10*L;
A=4;
rL=N/A;
lL=N-N/A+1;
unitPar1= sqrt([ones(rL,1); (1+cos((0:N)'*pi/N))/2; zeros(N/A-1,1)]);
unitPar2=flipud(unitPar1);

x=1:1:size(unitPar1);
k1=randomField( law, correl, L, mu, s, Nmc, x);
k2=randomField( law, correl, L, mu, s, Nmc, x);

vk3 = 0;
mk3 = 0;
for n=1:nSamples
    k1=randomField( law, correl, L, mu, s, Nmc, x);
    k2=randomField( law, correl, L, mu, s, Nmc, x);
    k3=k1.*unitPar1+k2.*unitPar2;
    vk3 = vk3+var(k3);
    mk3 = mk3+mean(k3);

end

vk3 = vk3/nSamples;
mk3 = mk3/nSamples;
%% Plot

rL=N/A;
lL=N/A+N+1;
fSize = 10;
fSize2 = 20;
lWeight = 2;

F1=figure; clf
set(F1,'defaulttextinterpreter','latex')
subplot(3,4,1:2)
plot(k1,'b', 'LineWidth',lWeight);
t=title('$$S1$$');
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
h=legend('S1', 'Location', 'southeast');
set(h,'interpreter','latex','fontsize',fSize)
set(t,'interpreter','latex','fontsize',fSize2)

subplot(3,4,3:4)
plot(k2,'r', 'LineWidth',lWeight);
t=title('$$S2$$');
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
h=legend('S2', 'Location', 'southeast');
set(h,'interpreter','latex','fontsize',fSize)
set(t,'interpreter','latex','fontsize',fSize2)

subplot(3,4,5:6)
hold on
plot(unitPar1,'m', 'LineWidth',lWeight);
t=title('$$S1_{modif} = S1 \sqrt{\Psi_1}$$');
plot(k1,'b--');
plot(unitPar1.*k1,'b', 'LineWidth',lWeight);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
h=legend('$$\sqrt{\Psi_1}$$','S1', '$$S1_{modif}$$', 'Location', 'southeast');
set(h,'interpreter','latex','fontsize',fSize)
set(t,'interpreter','latex','fontsize',fSize2)
hold off

subplot(3,4,7:8)
hold on
plot(unitPar2,'g', 'LineWidth',lWeight);
t=title('$$S2_{modif} = S2 \sqrt{\Psi_2}$$');
plot(k2,'r--');
plot(unitPar2.*k2, 'r', 'LineWidth',lWeight);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
h=legend('$$\sqrt{\Psi_2}$$','$$S2$$', '$$S2_{modif}$$', 'Location', 'southeast');
set(h,'interpreter','latex','fontsize',fSize)
set(t,'interpreter','latex','fontsize',fSize2)
hold off

subplot(3,4,10:11)
hold on
plot(unitPar1.*k1,'b--');
t=title('$$Resultant~Sample$$')
plot(unitPar2.*k2,'r--');
plot(unitPar1.*k1+unitPar2.*k2,'k', 'LineWidth',lWeight);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
h=legend('$$S1_{modif}$$','$$S2_{modif}$$','$$S1_{modif}+S2_{modif}$$', 'Location', 'southeast');
set(h,'interpreter','latex','fontsize',fSize)
set(t,'interpreter','latex','fontsize',fSize2)
hold off