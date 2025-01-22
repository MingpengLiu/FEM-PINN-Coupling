function [xdis, ydis]=dis_cal(Node_information,node_cycle,dis)
sintheta=Node_information(node_cycle,2)/0.2;
costheta=Node_information(node_cycle,3)/0.2;
xdis=sintheta*dis;
ydis=costheta*dis;
end