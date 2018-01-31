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
%nTimeSteps = 868;

%figure(2)
%plot(t_vec,P);

%Sensor coordinates (relative to source)
R_mat=[1, 0, 0; ...
       2, 0, 0; ...
       3, 0, 0];
%R_mat=[3, 0, 0];
nCoord=size(R_mat,1);

%Calculate displacement evolution for each R
distance = zeros(nCoord,1);
u_Rx=zeros(nCoord,nTimeSteps);
u_signal=zeros(nCoord,nTimeSteps);
PSI=zeros(nCoord,nTimeSteps);
X=zeros(nCoord,nTimeSteps);
Q=zeros(nCoord,nTimeSteps); 

for coord=1:nCoord

    R_vec = R_mat(coord,:);
    R = norm(R_vec);
    distance(coord) = R;
    t_p=R/alpha
    t_s=R/beta
    
 %TEST   
%     P = rickerSEM(t_vec,2,0.5,1);
%     mask = (t_vec>=t_p) .* (t_vec<=t_s);
%     P_trunc=P.*mask;
%     t_vec_trunc=t_vec.*mask;
%     %VERIFIER
%     Q(coord,:)=beta^2/R^2*ifft(conj(fft(t_vec_trunc)).*(fft(P_trunc)));
%     P_tp = rickerSEM(t_vec-t_p,2,0.5,1);
%     P_ts = rickerSEM(t_vec-t_s,2,0.5,1);
%     PSI(coord,:)= P_ts - Q(coord,:);
%     X(coord,:)=(((beta/alpha)^2)*P_tp)-P_ts+3*Q(coord,:);
%     u_Rx(coord,:) = (PSI(coord,:) + X(coord,:))/(4*pi*mu*R);
%     u_signal(coord,:) = u_Rx(coord,:);
    
    
  %OLD   
    t_s_pos=sum(t_vec<=t_s);
    P = rickerSEM(t_vec,2,0.5,1);
    %P = circshift(P,-t_s_pos);
    %P_tp = rickerSEM(t_vec-t_p,2,0.5,1);
    %P_ts = rickerSEM(t_vec-t_s,2,0.5,1);
    %P=P_ts;
    dirac_p = 1*(t_vec==t_p);
    dirac_s=1*(t_vec==t_s);
    hside_p=1/2*(t_vec==t_p)+(t_vec>t_p);
    hside_s=1/2*(t_vec==t_s)+(t_vec>t_s);
    
    PSI(coord,:)=dirac_s - t_vec.*(((beta/alpha)^2.*hside_p/t_p^2)-hside_s/t_s^2);
    X(coord,:)=(((beta/alpha)^2)*dirac_p-dirac_s)+3*t_vec.*(((beta/alpha)^2).*hside_p/t_p^2 -hside_s/t_s^2);
    u_Rx(coord,:) = (PSI(coord,:) + X(coord,:))/(4*pi*mu*R);
    conv = ifft(conj(fft(P)).*(fft(u_Rx(coord,:))));
    u_signal(coord,:) = conv;
    %u_signal(coord,:) = ifft(conj(fft(P)).*(fft(u_Rx(coord,:))));
    %u_signal(coord,:) = circshift(conv,+t_s_pos);
    
end

figure(3)
subplot(2,3,1);
hold on
for coord=1:nCoord
    plot(t_vec,PSI(coord,:))
    label{coord} = ['R = ',num2str(distance(coord))]; 
end
legend(label,'FontSize',20);
title('PSI','FontSize',20);
hold off

subplot(2,3,2);
hold on
for coord=1:nCoord
    plot(t_vec,X(coord,:))
    label{coord} = ['R = ',num2str(distance(coord))]; 
end
legend(label,'FontSize',20);
title('X','FontSize',20);
hold off

subplot(2,3,3);
hold on
for coord=1:nCoord
    plot(t_vec,u_Rx(coord,:))
    label{coord} = ['R = ',num2str(distance(coord))]; 
end
legend(label,'FontSize',20);
title('u_Rx','FontSize',20);
hold off

subplot(2,3,4);
hold on
for coord=1:1
    plot(t_vec,P)
    label{coord} = ['Input Signal']; 
end
legend(label,'FontSize',20);
title('Signal','FontSize',20);
hold off

subplot(2,3,5);
hold on
for coord=1:1
    plot(t_vec,fft(P))
    label{coord} = ['FFT Input Signal']; 
end
legend(label,'FontSize',20);
title('FFT Signal','FontSize',20);
hold off

subplot(2,3,6);
hold on
for coord=1:nCoord
    plot(t_vec,u_signal(coord,:))
    label{coord} = ['R = ',num2str(distance(coord))]; 
end
legend(label,'FontSize',20);
title('u_{signal}','FontSize',20);
hold off

% figure (1)
% %X
% %subplot(1,3,1);
% hold on
% count=0;
% for i = 1:size(u_signal,1)
%     count = count + 1;
%     plot(t_vec, u_signal(i,:));
%     label{count} = ['Analytic R = ',num2str(distance(i))];
% end
% title('Displacement X','FontSize',20);
% xlabel('Time [s]','FontSize',20);
% ylabel('Displacement [m]','FontSize',20);
% set(gca,'fontsize',20)
% hold off


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