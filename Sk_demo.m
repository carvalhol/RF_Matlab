
k = 0.1*(0:100);
x = 1;
Sk = exp(-k.^2*x);
cosk = cos(k.*x);
cosk2 = cos(k.*2);
lWeight = 5;
axisSize = 20;

hold on
plot(k, Sk,'b', 'LineWidth',lWeight);
plot(k, cosk,'m--', 'LineWidth',4);
%plot(k, cosk2,'g--', 'LineWidth',4);
box on
set(gca,'fontsize', axisSize);
hold off