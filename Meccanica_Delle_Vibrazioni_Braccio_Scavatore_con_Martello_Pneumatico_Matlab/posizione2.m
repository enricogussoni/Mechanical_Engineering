function f = posizione2(x)

global l
global theta2_tmp theta10

k = x(1);
tau = x(2);

m = 1;
eps = 0;

f=[m*cos(eps)+l*cos(theta10+theta2_tmp)-k*cos(tau);
    m*sin(eps)+l*sin(theta10+theta2_tmp)-k*sin(tau)];


