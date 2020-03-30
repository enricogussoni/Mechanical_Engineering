

%************** Ese 2: equazione del calore *****************************

clear all; close all
N = 100;
n = 1:N;

%Dato iniziale: 

fi = @(x) sin(x) + 0.5*sin(2*x) + 0.2*sin(5*x);

% calcolo coefficienti bk tramite dst
xs = pi*n/(N+1);
bk = 2/(N+1)*dst(fi(xs));

%Ricostruire la soluzione in spazio e tempo

% PUNTO 1: Soluzione ad un istante temporale fissato T=0,0.5,1

T=0; %Ritroviamo la condizione iniziale?
xdis = linspace(0,pi,1000);

u=0*xdis;
for i = n
    u = u + bk(i).*exp(-i.^2*T).*sin(i*xdis);
end

figure(1)
plot(xdis,fi(xdis),'ko', xdis, u, 'r'); 

%
T=0.5;
xdis = linspace(0,pi,1000);

u=0;
for i = n
    u = u + bk(i).*exp(-i.^2*T).*sin(i*xdis);
end

hold on
plot(xdis, u, 'y'); 

%
T=1;
xdis = linspace(0,pi,1000);

u=0;
for i = n
    u = u + bk(i).*exp(-i.^2*T).*sin(i*xdis);
end

hold on
plot(xdis, u, 'g'); 
legend('Condizione iniziale','T=0','T=0.5','T=1');

%%
%PUNTO 2: Soluzione al variare del tempo
close all
figure(2)

T=[0:0.05:5];

xdis = linspace(0,pi,1000);

for t = T
     u=0*xdis;
     for i = n
         u = u + bk(i).*exp(-i.^2*t).*sin(i*xdis);
     end
     plot(xdis,u,'b'); hold on
     drawnow
  
 end


%%
%PUNTO 3 : cambiamo la soluzione iniziale con una discontinua 

fi= @(x) (x>pi/4).*(x<3/4*pi);

% calcolo coefficienti bk tramite dst
xs = pi*n/(N+1);
bk = 2/(N+1)*dst(fi(xs));

figure(3)
T = linspace(0,1,100);
xdis = linspace(0,pi,1000);

for t = T
     u=0*xdis;
     for i = n
         u = u + bk(i).*exp(-i.^2*t).*sin(i*xdis);
     end
     plot(xdis,u,'b'); hold on
     drawnow
 end
