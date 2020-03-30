function [k,dk,ddk,tau,dtau,ddtau]=cinematica2(theta2)
global  l 
global theta2_tmp theta10
global k_old  tau_old
k_old=1.6;
tau_old=pi/6;
theta2_tmp=theta2;
l=0.5;
options=optimset('TolFun',1e-6,'Display','off');
x02 = [k_old tau_old];
x2 = fsolve(@posizione2,x02,options);

k = x2(1);
tau=x2(2);
Jaco=[cos(tau),-k*sin(tau);sin(tau),k*cos(tau)]\[-l*sin(theta10+theta2);l*cos(theta10+theta2)];
dk=Jaco(1);
dtau=Jaco(2);

Hes=[cos(tau) -k*sin(tau);sin(tau) k*cos(tau)]\[l*cos(theta10+theta2)-2*dk*dtau*sin(tau)-k*dtau^2*cos(tau);l*sin(theta10+theta2)+2*dk*dtau*cos(tau)-k*dtau^2*sin(tau)];
ddk=Hes(1);
ddtau=Hes(2);
k_old=k;
tau_old=tau;