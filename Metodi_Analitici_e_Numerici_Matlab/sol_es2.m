
clc; clear all; close all;
L1 = -100;
L2 = -1;

A = [0 1; -(L1*L2) L1+L2];
f = @(t,y) A*y;
y0 = [1;1];

[t,u] = eulero_avanti_sys(f,5,y0,0.0001);
figure; plot(t,u(1,:),'bo-',t,u(2,:),'ro-')
legend('y_1','y_2')

%% la condizione di assoluta stabilita' e' h<1/50
close all
[t1,u1] = eulero_avanti_sys(f,5,y0,0.019); 
figure; plot(t1,u1(1,:),'bo-',t1,u1(2,:),'ro-')
title('h = 0.019')
legend('u_1','u_2')

[t2,u2] = eulero_avanti_sys(f,5,y0,0.021); 
figure; plot(t2,u2(1,:),'bo-',t2,u2(2,:),'ro-')
title('h = 0.021')
legend('u_1','u_2')

[t2,u2] = ode45(f,[0 5],y0');
figure; plot(t2,u2(:,1),'bo-',t2,u2(:,2),'ro-'); title('ODE45');
figure; plot(t2(1:end-1),diff(t2),'.-'); title('ODE45');
length(t2)

[t3,u3] = ode15s(f,[0 5],y0');
length(t3)
figure; plot(t3,u3(:,1),'bo-',t3,u3(:,2),'ro-'); title('ODE15s');
figure; plot(t3(1:end-1),diff(t3),'.-'); title('ODE15s');
