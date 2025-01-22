figure(1)
plot(plotstrain(:,2),q,'r')
title('Deviatoric stress kPa');
hold off

figure(2)
plot(p,q,'r',[0 100],1.2*[0 100],'k')
title('Deviatoric stress kPa');
hold off


figure(3)
plot(plotstrain(:,2),epd,'k--')
legend('plastic deviatoric strain')
title('Strain');
hold off

figure(4)
plot(plotstrain(:,2),plote,'r')
title('Void ratio');
hold off