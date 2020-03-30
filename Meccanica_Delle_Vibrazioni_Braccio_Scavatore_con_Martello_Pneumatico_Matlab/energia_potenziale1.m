function [dy1,ddy1,dy2,ddy2,dy3,ddy3,DL,dDL,ddDL] = energia_potenziale1(theta)

global a s p x
global g0

dy1 = a*cos(theta);
ddy1 = -a*sin(theta); 

dy2=(a+s)*cos(theta)+p*cos(theta-(10*pi)/24);
ddy2= -(a+s)*sin(theta)-p*sin(theta-(10*pi)/24);

dy3=(a+s)*cos(theta)+(p+x)*cos(theta-(10*pi)/24);
ddy3= -(a+s)*sin(theta)-(p+x)*sin(theta-(10*pi)/24);

[g,dg,ddg] = cinematica1(theta);

DL = g-g0;
dDL = dg;
ddDL = ddg;
