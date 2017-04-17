close all
clear all

set(0,'defaulttextinterpreter','latex')

L = 2000;
partition = ones(L);

ovlp = 300;
bord = L - ovlp + 1;

phi = ones(ovlp,L);
for i =1:ovlp
    phi(i,:)= phi(i,:)*(i-1)/ovlp;
end


partition(1:ovlp,:) = partition(1:ovlp, :).*phi(:,:);
partition(bord:L,:) = partition(bord:L, :).*flip(phi(:,:),1);
partition(:,1:ovlp) = partition(:,1:ovlp).*transpose(phi(:,:));
partition(:,bord:L) = partition(:,bord:L).*transpose(flip(phi(:,:),1));


figure(1)
surf(partition,'EdgeColor','none','LineStyle','none')
set(gca,'FontSize',20)
axis off
box on
%colormap(gray)
colormap(jet)