clear all
close all
clc


%nlhet_setup_postprocess;
%nlhet_load_traces;
load('/Users/carvalhol/Desktop/coherenceAng/case_1_and_case_2.mat')

%% Set by user
case_Nb = 1;
t_max = 10; %Ne pas changer (voir dans les snapshots quand ?a va a zero
t=[0.:0.001:t_max];
%Temps pour la phase forte
t_scope_min = 2.;
t_scope_max = 5.;
cp = 3; %captors Step. Each step is 10 meters
capt_Max=534; %Last sensor to be used in analisys
capt_Nb=[  1:41-cp   42:82-cp   83:123-cp 124:164-cp...
         165:205-cp 206:246-cp 247:287-cp 288:328-cp...
         329:369-cp 370:410-cp 411:451-cp 452:492-cp...
         493:533-cp]';%sensors divided in intervals according to distribution

tt=sem{case_Nb}.Time(:,1)';
xx=sem{case_Nb}.Displ.x';
xy=sem{case_Nb}.Displ.y';
xz=sem{case_Nb}.Displ.z';
dt=t(2)-t(1);
nt_1=int16(t_scope_min/dt);
nt_2=int16(t_scope_max/dt);

%-------------- discretisation
M  = 5; %Abrahamson
OM = pi/dt;  % frequence de coupure
Nf = 1 + (t_scope_max - t_scope_min)/dt; %nombre de discr?tisation pour la phase forte
dw = 2*OM/Nf ; % step de frequence
df = dw/2./pi ; %Hz

w=0:dw:(Nf-1)*dw;
fr=w/2./pi;
f_f=fr(1:int16(Nf/2)-2*M);
f=[0.:0.1:20];

ddw=dw;

% %ANG START
% cd '.';
% 
% data=dir;
% da=importdata(data(1+3).name);
% tt=da(:,1);
% xx(1,:)=da(:,2);
% xy(1,:)=da(:,3);
% xz(1,:)=da(:,4);
% 
% t=[0.:0.0005:7];
% t_scope_min = 1.2;
% t_scope_max = 5.;
% capt_Max = 40;
% 
% capt_Nb=[1:10 12:20 22:30 32:39]';
% cp = 1;
% 
% %ANG END

%Reading Files

%Processing Data
for c=1:capt_Max
    
%     %ANG START
%     data=dir;
%     da=importdata(data(c+3).name);
%     tt=da(:,1);
%     xx(c,:)=da(:,2);
%     xy(c,:)=da(:,3);
%     xz(c,:)=da(:,4);
%     dt=t(2)-t(1);
%     nt_1=int16(t_scope_min/dt);
%     nt_2=int16(t_scope_max/dt);
% 
%     OM=pi/dt;  % frequence de coupure
%     Nf=length(pfx(c,:)); %nombre de discr?tisation pour la phase forte
%     dw=2*OM/Nf ; % step de frequence
%     df=dw/2./pi ; %Hz
%     
%     w=0:dw:(Nf-1)*dw;
%     fr=w/2./pi;
%     f_f=fr(1:int16(Nf/2)-2*M);
%     f=[0.:0.1:20];
%     
%     ddw=dw;
%     %ANG END
    
    xa(c,:)=interp1(tt,xx(c,:),t);
    ya(c,:)=interp1(tt,xy(c,:),t);
    za(c,:)=interp1(tt,xz(c,:),t);
    
    %Phase forte
    pfx(c,:)=xa(c,[nt_1:nt_2]);
    pfy(c,:)=ya(c,[nt_1:nt_2]);
    pfz(c,:)=za(c,[nt_1:nt_2]);
   
    pfx(c,:)=cosine_taper(pfx(c,:),0.05);
    pfy(c,:)=cosine_taper(pfy(c,:),0.05);
    pfz(c,:)=cosine_taper(pfz(c,:),0.05);
         
    Ux(c,:)=fft(pfx(c,:));
    Uy(c,:)=fft(pfy(c,:));
    Uz(c,:)=fft(pfz(c,:));
    
    j_inf = M+1;
    j_sup = int16(Nf/2)-M;
    j_range = j_sup - j_inf + 1;
    
    if(c==1)
        SMux = zeros(capt_Max,j_range);
        SMuy = zeros(capt_Max,j_range);
        SMuz = zeros(capt_Max,j_range);
    end
    
    j_c = 0; %j counter
    
    for j=j_inf:j_sup
        
        j_c = j_c + 1;
        
        for m=-M:M
            
            W=(0.54 -0.46*cos((pi*((m*ddw)+M))/M))/(1.08*M);
            
            SMux(c,j_c) = SMux(c,j_c)+W*(conj(Ux(c,j+m)).*(Ux(c,j+m)));
            SMuy(c,j_c) = SMuy(c,j_c)+W*(conj(Uy(c,j+m)).*(Uy(c,j+m)));
            SMuz(c,j_c) = SMuz(c,j_c)+W*(conj(Uz(c,j+m)).*(Uz(c,j+m)));
            
        end
    end
    
    %keyboard
    
end

p = [capt_Nb capt_Nb+cp]; %captors pairs to be compared

coh10x1=zeros(size(p,1), numel(f));
coh10y1=coh10x1;
coh10z1=coh10x1;
SMCx10 =zeros(size(p,1), j_range);
SMCy10 =SMCx10;
SMCz10 =SMCx10;

for k=1:size(p,1)
    
    j_c = 0;
    
    for j=j_inf:j_sup
        
        j_c = j_c + 1;
        
        for m=-M:M
            
            W=(0.54 -0.46*cos((pi*((m*ddw)+M))/M))/(1.08*M);
            
            SMCx10(k,j_c) = SMCx10(k,j_c) + W*(conj(Ux(p(k,1),j+m)).*(Ux(p(k,2),j+m)));
            SMCy10(k,j_c) = SMCy10(k,j_c) + W*(conj(Uy(p(k,1),j+m)).*(Uy(p(k,2),j+m)));
            SMCz10(k,j_c) = SMCz10(k,j_c) + W*(conj(Uz(p(k,1),j+m)).*(Uz(p(k,2),j+m)));
            
        end
                
    end
    
    coh10x1(k,:)=interp1(f_f,(SMCx10(k,:)./(sqrt(SMux(p(k,1),:)).*sqrt(SMux(p(k,2),:)))),f);
    coh10y1(k,:)=interp1(f_f,(SMCy10(k,:)./(sqrt(SMuy(p(k,1),:)).*sqrt(SMuy(p(k,2),:)))),f);
    coh10z1(k,:)=interp1(f_f,(SMCz10(k,:)./(sqrt(SMuz(p(k,1),:)).*sqrt(SMuz(p(k,2),:)))),f);
    
end

% figure(1)
% hold on
% for k=1:size(p,1)
%     plot(f, abs(coh10x1(k,:)));
% end
% hold off

coherency10x=real(tanh(mean(atanh(coh10x1),1)));
coherency10y=real(tanh(mean(atanh(coh10y1),1)));
coherency10z=real(tanh(mean(atanh(coh10z1),1)));

figure(2)
plot(f,coherency10x,f,coherency10y)
ylim([0,1]);