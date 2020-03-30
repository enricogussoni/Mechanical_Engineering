function [t,x,xp]=rk4(statofun,tf,dt,x0)
% Integrazione con il metodo di Runge-Kutta del 4. Ordine a passo fisso


t=[0:dt:tf];
np=length(t);


kut=[dt/6;2*dt/6;2*dt/6;dt/6];

x(:,1)=x0;
xp(:,1)=feval(statofun,0,x0);

att= waitbar(0,'Attendere, integrazione numerica in corso...');

for i=2:np
    
    tempx(:,1)=feval(statofun,t(i-1),x(:,i-1));
    tempx(:,2)=feval(statofun,t(i-1)+dt/2,x(:,i-1)+dt/2*tempx(:,1));
    tempx(:,3)=feval(statofun,t(i-1)+dt/2,x(:,i-1)+dt/2*tempx(:,2));
    tempx(:,4)=feval(statofun,t(i-1)+dt,x(:,i-1)+dt*tempx(:,3));
    
    x(:,i)=x(:,i-1)+tempx*kut;
    xp(:,i)=feval(statofun,t(i),x(:,i));
    
    waitbar(i/np)
    
end

close (att);