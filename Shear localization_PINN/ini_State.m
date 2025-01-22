function [ini_stress,ini_e0] = ini_State(Element_information,Node_information,e01,e02,e03,p01,p02,p03)
h1=0.12;h2=0.04;

nEle = length(Element_information);
ini_stress = zeros(nEle,6);
ini_e0 = zeros(nEle,1);
for i=1:nEle
    node = Element_information(i,2:4);
    nx0_e = Node_information(node,2:3);
    if all(nx0_e(:,2)>h1)  
        ini_stress(i,:)=[-p01 -p01 -p01 0 0 0];
        ini_e0(i,:)=e01;
    elseif all(nx0_e(:,2)<h1 & nx0_e(:,2)>h2) 
        ini_stress(i,:)=[-p02 -p02 -p02 0 0 0];
        ini_e0(i,:)=e02;
    else 
        ini_stress(i,:)=[-p03 -p03 -p03 0 0 0];
        ini_e0(i,:)=e03;
    end
end