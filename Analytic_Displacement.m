close all
clear all
clc

%% Analitical Displacement
% Eduardo Kausel page 48
% Three-dimensional problems in full,homogeneous spaces
% 4.2 Point load (Stokes problem)

alpha = 5; %Vp
beta = alpha/sqrt(2); %Vs
rho = 2800; %Density
mu = rho*beta^2; %G

%T = 1.5; %period
T = 5; %period
dT = 0.001;
nTimeSteps = 1+ceil(T/dT);
t_vec=linspace(0,T,nTimeSteps);
Signal = rickerSEM(t_vec,2,0.5,1);
%nTimeSteps = 868;

figure(2)
plot(t_vec,Signal);

%Sensor coordinates (relative to source)
R_mat=[1, 0, 0; ...
       2, 0, 0; ...
       3, 0, 0];
%R_mat=[3, 0, 0];
distance = zeros(size(R_mat,1),1);

%Calculate displacement evolution for each R

u_signal=zeros(size(R_mat,1),nTimeSteps);

for coord=1:size(R_mat,1)

    R_vec = R_mat(coord,:);
    R = norm(R_vec);
    distance(coord) = R;
    t_p=R/alpha;
    t_s=R/beta;


    %Time evolution
    u_Rx=zeros(1,nTimeSteps);
    count=0;

    for count=1:nTimeSteps
        t_now=t_vec(count);
        dirac_p=1*(t_now==t_p);
        dirac_s=1*(t_now==t_s);
        hside_p=1/2*(t_now==t_p)+(t_now>t_p);
        hside_s=1/2*(t_now==t_s)+(t_now>t_s);

        PSI=dirac_s - t_now*(((beta/alpha)^2*hside_p/t_p^2)-hside_s/t_s^2);
        X=(beta/alpha)^2*dirac_p-dirac_s+3*t_now*(((beta/alpha)^2)*hside_p/t_p^2 -hside_s/t_s^2);

        u_Rx(count) = (PSI + X)/(4*pi*mu*R);

    end

    %Signal = ricker(t_vec,2,0.5);
    u_signal(coord,:) = ifft(conj(fft(Signal)).*(fft(u_Rx)));
end

figure (1)
%X
%subplot(1,3,1);
hold on
count=0;
for i = 1:size(u_signal,1)
    count = count + 1;
    plot(t_vec, u_signal(i,:));
    label{count} = ['Analytic R = ',num2str(distance(i))];
end
title('Displacement X','FontSize',20);
legend(label,'FontSize',20);
xlabel('Time [s]','FontSize',20);
ylabel('Displacement [m]','FontSize',20);
set(gca,'fontsize',20)
hold off


% %% 
% % Eduardo Kausel page 48
% % Three-dimensional problems in full,homogeneous spaces
% % 4.2Point load (Stokes problem)
% 
% alpha = 5; %Vp
% beta = alpha/sqrt(2); %Vs
% rho = 2800; %Density
% mu = rho*beta^2; %G
% 
% T = 1.5; %period
% nTimeSteps = 868;
% 
% %Sensor coordinates (relative to source)
% R_mat=[1, 0, 0; ...
%        2, 0, 0; ...
%        3, 0, 0];
% 
% %Calculate displacement evolution for each R
% figure(1)
% hold on
% for coord=1:size(R_mat,1)
% 
%     R_vec = R_mat(coord,:);
%     R = norm(R_vec);
%     t_p=R/alpha;
%     t_s=R/beta;
% 
% 
%     %Time evolution
%     u_Rx=zeros(1,nTimeSteps);
%     count=0;
%     t_vec=linspace(0,T,nTimeSteps);
% 
%     for count=1:nTimeSteps
%         t=t_vec(count);
%         dirac_p=1*(t==t_p);
%         dirac_s=1*(t==t_s);
%         hside_p=1/2*(t==t_p)+(t>t_p);
%         hside_s=1/2*(t==t_s)+(t>t_s);
% 
%         PSI=dirac_s - t*(((beta/alpha)^2*hside_p/t_p^2)-hside_s/t_s^2);
%         X=(beta/alpha)^2*dirac_p-dirac_s+3*t*(((beta/alpha)^2)*hside_p/t_p^2 -hside_s/t_s^2);
% 
%         u_Rx(count) = (PSI + X)/(4*pi*mu*R);
% 
%     end
% 
%     Signal = ricker(t_vec,2,0.5);
%     u_signal = ifft(fft(Signal).*fft(u_Rx));
%     plot(t_vec, u_signal);
% end
% hold off