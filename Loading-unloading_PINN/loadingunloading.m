function [disp1, disp2] = loadingunloading(strainpoints,numberpoints)

disp1=[];
disp2=[];
% first loading
for i=1:numberpoints(1)
    disp1(i) = strainpoints(1)/numberpoints(1)*i;
    disp2(i) = disp1(i);
end

% unloading
for i=1:numberpoints(2)
    disp1(end+1) = strainpoints(1)-abs(strainpoints(2)-strainpoints(1))/numberpoints(2)*i;
    disp2(end+1) = disp1(end);
end

% re-loading
for i =1:numberpoints(3)
    disp1(end+1) = strainpoints(2)+(strainpoints(3)-strainpoints(2))/numberpoints(3)*i;
    disp2(end+1) = disp1(end);
end

end