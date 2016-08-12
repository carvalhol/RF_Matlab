clc
close all
clear all

lWeight = 2;
axisSize = 2;
xStep = 0.1;

x = xStep*(0:100);
l_c = 1.0;
testN = 0;
L_norm = 2;

corrType =['E','P','G']; %E = Exponential, P = Powerlaw, G = Gaussian
phyType =['T','C']; %T = triangular, C = cossinus
alpha = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

results = zeros(numel(corrType)*numel(phyType), numel(alpha));

for c = 1:numel(corrType)
    
    for p = 1:numel(phyType)
        
        for a = 1:numel(alpha)

           testN = testN + 1;
           fprintf('---\n');
           fprintf('Test %d\n', testN);
           fprintf('alpha = %f\n', alpha(a));

           switch(corrType(c))
               case 'E'
                  fprintf('corrType: Exponential\n' );
                  R = exp(-abs(x)/l_c); %Exponential
                  if(a==1)
                      R_E = R;
                  end
               case 'P'
                  fprintf('corrType: Power Law\n' );
                  R = (1 + abs(x.^2)/l_c).^(-2); %Power Law
                  if(a==1)
                      R_P = R;
                  end
               case 'G' 
                  fprintf('corrType: Gaussian\n' );
                  R = exp(-abs(x.^2)/l_c^2); %Gaussian
                  if(a==1)
                      R_G = R;
                  end
               otherwise
                  fprintf('Invalid corrType\n' );
           end

           switch(phyType(p))
               case 'T'
                  fprintf('phyType: Triangular\n' );
                  phi = 1 - x/alpha(a)/l_c; phi(phi<0) = 0; %Triangular
               case 'C'
                  fprintf('phyType: Cosinus\n' );
                  phi = 1/2+(1/2*(cos(x*pi/alpha(a)/l_c))); phi(x/alpha(a)/l_c>1) = 0; %Cosinus
               otherwise
                  fprintf('Invalid phyType\n' );
           end

        %R = exp(-abs(x.^2)/l_c^2); %Gaussian
        %R = exp(-abs(x)/l_c); %Exponential
        %R = (1 + abs(x.^2)/l_c)^(-2); %Power Law


        %phi = 1 - x/alpha/l_c; phi(phi<0) = 0; %Triangular
        %phi = 1/2+(1/2*(cos(x*pi/alpha/l_c))); phi(x/alpha/l_c>1) = 0; %Cosinus

            R_overlap = sqrt(phi).*R;

            Dist_L2 = (trapz(x,(R-R_overlap).^L_norm))^(1/L_norm);
            %Dist_L2 = abs((trapz(x,(R).^L_norm))^(1/L_norm)-(trapz(x,(R_overlap).^L_norm))^(1/L_norm));
            Dist_L2_norm = Dist_L2/(trapz(x,R.^L_norm)^(1/L_norm));
            
            fprintf('Dist_L2 %f\n', Dist_L2);
            
            %results((c-1)*numel(phyType)+p,a) = Dist_L2;
            results((c-1)*numel(phyType)+p,a) = Dist_L2_norm;
            
            
        end
    end
end

figure(1)
l_weigth = 4;
hold all
plot(R_E, 'LineWidth',l_weigth)
plot(R_P, 'LineWidth',l_weigth)
plot(R_G, 'LineWidth',l_weigth)
Legend = cell(3,1);
Legend(1) = {'Exponential'};
Legend(2) = {'Power Law'};
Legend(3) = {'Gaussian'};
xlim([0,60]);
set(gca,'xtick',[],'ytick',[0, 1], 'LineWidth', l_weigth, 'FontSize', 40);

set(1, 'Position', [0, 0, 550, 500]);
rule_fig(1);
leg=legend(Legend,'Location','northeast','FontSize',40);
set(leg,'Interpreter','latex');
saveas(1,'Correlation_Functions','epsc');

% R_overlap = zeros(1, numel(R));
% count_dist = zeros(1, numel(R));
% for i = 1:numel(x)
%     for j = 1:numel(x)
%         dist = abs(j-i)+1;
%         R_overlap(dist) = R_overlap(dist)+(phi(i)*phi(j));
%         count_dist(dist) = count_dist(dist)+1;
%     end
% end
% R_overlap = sqrt(R_overlap./count_dist).*R;

% %% Plot
% hold on
% plot(x, R, 'LineWidth',lWeight);
% plot(x, phi, 'LineWidth',lWeight);
% plot(x, R_overlap, 'LineWidth',lWeight);
% %plot(k, cosk,'m--', 'LineWidth',4);
% %plot(k, cosk2,'g--', 'LineWidth',4);
% box on
% set(gca,'fontsize', axisSize);
% hold off