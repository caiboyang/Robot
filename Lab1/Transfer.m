function [T04] = Transfer(q)
T01=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 q(1); 0 0 0 1];
T12=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*[1 0 0 0; 0 1 0 0; 0 0 1 q(2); 0 0 0 1];
T23=[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]*[cos(q(3)) -sin(q(3)) 0 0; sin(q(3)) cos(q(3)) 0 0; 0 0 1 0; 0 0 0 1];
T34=[1 0 0 1; 0 0 -1 0; 0 1 0 0; 0 0 0 1]*[cos(q(4)) -sin(q(4)) 0 0; sin(q(4)) cos(q(4)) 0 0; 0 0 1 0; 0 0 0 1];
T04=T01*T12*T23*T34;