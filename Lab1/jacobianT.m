function [T24,y] = jacobianT(q,p)
T01=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 q(1); 0 0 0 1];
T12=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 q(2); 0 0 0 1];
T23=[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]*[cos(q(3)) -sin(q(3)) 0 0; sin(q(3)) cos(q(3)) 0 0; 0 0 1 0; 0 0 0 1];
T34=[1 0 0 1; 0 0 -1 0; 0 1 0 0; 0 0 0 1]*[cos(q(4)) -sin(q(4)) 0 0; sin(q(4)) cos(q(4)) 0 0; 0 0 1 0; 0 0 0 1];
% T04=T01*T12*T23*T34;
T24=T23*T34;

% y=[ 0, 0, p4(3)*cos(q(3)) - p4(1)*cos(q(4))*sin(q(3)) + p4(2)*sin(q(3))*sin(q(4)), - p4(2)*cos(q(3))*cos(q(4)) - p4(1)*cos(q(3))*sin(q(4));
%  0, 0,                                                       0,                   p4(1)*cos(q(4)) - p4(2)*sin(q(4));
%  0, 0, p4(2)*cos(q(3))*sin(q(4)) - p4(1)*cos(q(3))*cos(q(4)) - p4(3)*sin(q(3)),   p4(2)*cos(q(4))*sin(q(3)) + p4(1)*sin(q(3))*sin(q(4));
%  0, 0,                                                       0,                                           0];
y= [0, 0, p(3)*cos(q(3)) - p(4)*sin(q(3)) - p(1)*cos(q(4))*sin(q(3)) + p(2)*sin(q(3))*sin(q(4)), - p(2)*cos(q(3))*cos(q(4)) - p(1)*cos(q(3))*sin(q(4));
 0, 0,                                                                 0,                   p(1)*cos(q(4)) - p(2)*sin(q(4));
 0, 0, p(2)*cos(q(3))*sin(q(4)) - p(3)*sin(q(3)) - p(1)*cos(q(3))*cos(q(4)) - p(4)*cos(q(3)),   p(2)*cos(q(4))*sin(q(3)) + p(1)*sin(q(3))*sin(q(4));
 0, 0,                                                                 0,                                         0];
end

