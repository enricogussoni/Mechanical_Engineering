clc; close all; clear all;
fun=@(t,u,v) -26*u - 2*v;
dfun = @(t,u,v) -2;
u0 = 1;
v0 = 1;

dt = 0.02;
T = 5;
[t1,s1,v1] = leapfrog(fun,dfun,T,u0,v0,dt);


plot(t1,s1,'bo-')
legend('leapfrog')

u_ex = @(t) exp(-t).*((2/5)*sin(5*t) + cos(5*t));
t_dis = linspace(0,T,1000);
hold on; plot(t_dis,u_ex(t_dis),'k')

%%
H = [0.2 0.1 0.05 0.025 0.0125 0.00625];
E1 = [];
for dt = H
    [t1,s1,v1] = leapfrog(fun,dfun,T,u0,v0,dt);
    E1 = [E1 max(abs(s1-u_ex(t1)))];
end
figure
loglog(H,E1,'bo-',H,H,'k',H,H.^2,'k--')
legend('Leapfrog', 'h', 'h^2')
