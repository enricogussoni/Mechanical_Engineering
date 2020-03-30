function [g,dg,ddg,delta,ddelta,dddelta] = cinematica1(theta)

global a 
global theta_tmp
global g_old delta_old 

theta_tmp = theta;

options = optimset('TolFun',1e-6,'Display','off');
x0 = [delta_old g_old];
x = fsolve(@posizione1,x0,options);

delta = x(1);
g=x(2);

J = [sin(delta),g*cos(delta);cos(delta),-g*sin(delta)]\[a*cos(theta);-a*sin(theta)];
ddelta =J(1);
dg =J(2); 

H = [sin(delta),g*cos(delta);cos(delta),-g*sin(delta)]\[a*sin(theta)+2*dg*ddelta*cos(delta)-ddelta^2*g*sin(delta);a*cos(theta)-2*dg*ddelta*sin(delta)-ddelta^2*g*cos(delta)];
dddelta = H(1);
ddg = H(2);

g_old = g;
delta_old = delta;