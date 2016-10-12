%Forward Kinematics
% unit for d1 and d2 are inch, unit for theta3 and theta4 are degree
% 0<d1<20
d1 = 10;
% 0<d2<11
d2 = 6;
theta3=30/180*pi;
theta4=30/180*pi;
% setup p4
p4=[0.817;1.366;16.683;1];

T01=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 d1; 0 0 0 1];
T12=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 d2; 0 0 0 1];
T23=[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]*[cos(theta3) -sin(theta3) 0 0; sin(theta3) cos(theta3) 0 0; 0 0 1 0; 0 0 0 1];
T34=[1 0 0 1; 0 0 -1 0; 0 1 0 0; 0 0 0 1]*[cos(theta4) -sin(theta4) 0 0; sin(theta4) cos(theta4) 0 0; 0 0 1 0; 0 0 0 1];

T04=T01*T12*T23*T34;
p0=T04*p4;



%reverse kinematics
xactual=[8.3628;1.5915;30.4356;1];
p4=[5;5;16;1];
q=[0;0;pi/4;pi/4];
[TT24,jacT]=jacobianT(q,p4);
x=TT24*p4; 
deltax = xactual-x;
a=sqrt((deltax(1)^2)+(deltax(2)^2));
while a>0.05
    deltaq=transpose(jacT)*deltax/5;
    q(3)=mod(deltaq(3)+q(3),2*pi);
    q(4)=mod(deltaq(4)+q(4),2*pi);
    [TT24,jacT]=jacobianT(q,p4);
    x=TT24*p4; 
    deltax = xactual-x;
    a=sqrt((deltax(1)^2)+(deltax(2)^2));
end
deltaz=xactual(3)-x(3);
    if (q(2)+deltaz)>11
    q(1)=q(1)+q(2)+deltaz-11;
    q(2)=11;
    elseif (q(2)+deltaz)<0
    q(1)=q(1)+(q(2)+deltaz);
    q(2)=0;
    else
    q(2)=q(2)+deltaz;
    end
T04=Transfer(q);
x=T04*p4;