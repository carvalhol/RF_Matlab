plot(1:10)

set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
set(gca,'Visible','off','plotboxaspectratio',[1 1 1]);
print -depsc test_print