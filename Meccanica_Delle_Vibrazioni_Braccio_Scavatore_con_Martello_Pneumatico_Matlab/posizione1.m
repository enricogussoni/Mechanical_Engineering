function f = posizione1(x)

global a 
global theta_tmp

g = x(2);
delta = x(1);

h = 1;
xi = -pi/4;

f = [h*cos(xi)+g*cos(delta)-a*cos(theta_tmp-pi/15); 
    h*sin(xi)+g*sin(delta)-a*sin(theta_tmp-pi/15)];
   

