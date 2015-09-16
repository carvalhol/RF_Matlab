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

F1=figure; clf
set(F1,'defaulttextinterpreter','latex')
subplot(3,4,1:2)
plot(k1,'b', 'LineWidth',2);
title('$$Delta$$','interpreter','latex')
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')

subplot(3,4,3:4)
plot(k2,'r', 'LineWidth',2);
title('Sample 2')
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')

subplot(3,4,5:6)
hold on
plot(unitPar1,'m', 'LineWidth',2);
plot(k1,'b--');
plot(unitPar1.*k1,'b', 'LineWidth',2);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
legend('\sqrt(2)Partition of Unity', 'Location', 'southeast')
hold off

subplot(3,4,7:8)
hold on
plot(unitPar2,'g', 'LineWidth',2);
plot(k2,'r--');
plot(unitPar2.*k2, 'r', 'LineWidth',2);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
legend('Partition of Unity', 'Location', 'southeast')
hold off

subplot(3,4,10:11)
hold on
plot(unitPar1.*k1,'b--');
title('Resultant Sample')
plot(unitPar2.*k2,'r--');
plot(unitPar1.*k1+unitPar2.*k2,'k', 'LineWidth',2);
ylim([-3,3]);
line([rL rL],ylim,'Color','k')
line([lL lL],ylim,'Color','k')
hold off