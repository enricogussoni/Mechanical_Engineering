function [dy4,ddy4,dy5,ddy5,dy6,ddy6,DL2,dDL2,ddDL2] = energia_potenziale2(theta1,theta2)

global a s p x 
global k0

dy4 = a*cos(theta1);   
ddy4 = -a*sin(theta1);

dy5 = (a+s)*cos(theta1)+ p*cos(theta1+theta2+7/6*pi);       
ddy5 = -(a+s)*sin(theta1)- p*sin(theta1+theta2+7/6*pi);  

dy6 = (a+s)*cos(theta1)+(p+x)*cos(theta1+theta2+7/6*pi);     
ddy6 = -(a+s)*sin(theta1)-(p+x)*sin(theta1+theta2+7/6*pi);

[k,dk,ddk] = cinematica2(theta2);

DL2 = k-k0;
dDL2 = dk;
ddDL2 = ddk;
