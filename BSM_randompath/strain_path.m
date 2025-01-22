function strainpath = strain_path(ns, nc, strainrange)
space = ns/(nc-1);
location=zeros(1,1);
for nn =1:nc-1
    location(end+1,1)=location(end,1)+space;
end

kesai11=zeros(1,1);
for i = 1:nc-1
   point = strainrange(1)+(strainrange(2)-strainrange(1))*rand(1);
   kesai11(end+1,1)=point;
end
kesai11=[location,kesai11];
kesai11 = interp1(kesai11(:,1), kesai11(:,2),1:ns, "cubic");
kesai11 = [0;kesai11'];

kesai22=zeros(1,1);
for i = 1:nc-1
   point = strainrange(1)+(strainrange(2)-strainrange(1))*rand(1);
   kesai22(end+1,1)=point;
end
kesai22=[location,kesai22];
kesai22 = (interp1(kesai22(:,1), kesai22(:,2),1:ns, "cubic"))';
kesai22 = [0;kesai22];

kesai12=zeros(1,1);
for i = 1:nc-1
   point = strainrange(1)+(strainrange(2)-strainrange(1))*rand(1);
   kesai12(end+1,1)=point;
end
kesai12=[location,kesai12];
kesai12 = (interp1(kesai12(:,1), kesai12(:,2),1:ns, "cubic"))';
kesai12 = [0;kesai12];

strainpath = [kesai11, kesai22, zeros(ns+1,1), kesai12, zeros(ns+1,1), zeros(ns+1,1)];
end