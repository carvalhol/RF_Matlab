clc
clear all
close all

xMax = 10;
deltaX_ref = 0.001;
x_ref = 0:deltaX_ref:xMax;
R_theoretical_1 = exp(-(x_ref).^2);
R_theoretical_2 = exp(-(x_ref));

deltaX_sample=1:-0.01:deltaX_ref;
epsilonD = zeros(numel(deltaX_sample),2);

for i=1:numel(deltaX_sample)
    x_sample = floor(x_ref/deltaX_sample(i))*deltaX_sample(i);
    
    %Gaussian
    R_sample = exp(-(x_sample).^2);
    bias = R_theoretical_1 - R_sample;
    epsilonD(i,1) = sqrt((sum(abs(bias).^2))*deltaX_ref);
    
    %Exponential
    R_sample = exp(-(x_sample));
    bias = R_theoretical_2 - R_sample;
    epsilonD(i,2) = sqrt((sum(abs(bias).^2))*deltaX_ref);
end

hold on
hold all

kMax = 2*pi./deltaX_sample;
deltaK = 2*pi/xMax;
n_K = kMax./deltaK;
plot(n_K, epsilonD(:,1),'k:','LineWidth', 2);
plot(n_K, epsilonD(:,2),'k--','LineWidth', 2);
legend({'Gaussian','Exponential'},'FontSize',20,'FontWeight','bold')
xlabel('N_k', 'FontSize', 20,'FontWeight','bold');
ylabel('Epsilon_D', 'FontSize', 20,'FontWeight','bold')
set(gca,'xscale','log', 'yscale','log', 'FontSize',15)
hold off