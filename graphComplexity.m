Vsize = 20;
y = zeros(Vsize, 1);
x = 1:Vsize;

for i = 1 : Vsize
    y(i) = theoretical_complexity([0.1, 0.1, 0.1], [i,i, i], 'gaussian', 10, 1.1);
end

plot(x, y, 'LineWidth', 3)

xlabel('(L/l_c)', 'FontSize', 20);
ylabel('Number of operations', 'FontSize', 20) 
set(gca,'FontSize',15)

grid('on')
box('on')