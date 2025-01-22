
figure(3);
clf;
subplot(2,2,1)
numstrain = 1:length(strainpath);
plot([numstrain',numstrain',numstrain',numstrain'],strainpath(:,1:4),"LineWidth",2);
ylabel("strain")

subplot(2,2,2)
numstress=1:length(save_stress);
plot([numstress',numstress',numstress',numstress'],save_stress(:,1:4), "LineWidth",2)
ylabel("stress")

subplot(2,2,3)
numstrainp=1:length(save_strain_p);
plot([numstrainp',numstrainp',numstrainp',numstrainp'],save_strain_p(:,1:4), "LineWidth",2)
ylabel("plastic strain")

subplot(2,2,4)
numvoid=1:length(save_void);
plot(numvoid',save_void(:), "LineWidth",2)
ylabel("void ratio")