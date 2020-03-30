close all; clear all; clc

GM = 398600; 
RT = 6360;
RL = 1737;
ang = [0:pi/100:2*pi];

twobody = @(t,y) [y(3); y(4); -GM*y(1)/(y(1)^2+y(2)^2)^(3/2); -GM*y(2)/(y(1)^2+y(2)^2)^(3/2)];

X0 = 405500;
V0 = 0.964;

day = 24*60*60;
[t45,u45] = ode45(twobody,[0 365*day], [X0 0 0 V0]);
fprintf('1) Numero di timestep: %d\n',length(t45));
figure(1); 
subplot(121); plot(u45(:,1),u45(:,2),'b-'); xlabel('x'); ylabel('y');
hold on; patch(RT*sin(ang),RT*cos(ang),'r'); 
subplot(122); plot(t45,sqrt(u45(:,1).^2+u45(:,2).^2)/X0,'b'); xlabel('t'); ylabel('|r|');
%% punto 2

opts = odeset('reltol',1e-6);
[t45,u45] = ode45(twobody,[0 365*day], [X0 0 0 V0], opts);
fprintf('2) Numero di timestep cambiando tolleranza: %d\n',length(t45));
figure(2); 
subplot(121); plot(u45(:,1),u45(:,2),'b'); xlabel('x'); ylabel('y');
hold on; patch(RT*sin(ang),RT*cos(ang),'r'); 
subplot(122); plot(t45,sqrt(u45(:,1).^2+u45(:,2).^2),'b'); xlabel('t'); ylabel('|r|');

%% punto 3
C = 1e-6;
twobody2 = @(t,y) [y(3); y(4); -GM*y(1)/(y(1)^2+y(2)^2)^(3/2)-C*y(3); -GM*y(2)/(y(1)^2+y(2)^2)^(3/2)-C*y(4)];

opts = odeset('reltol',1e-6);
[t45,u45] = ode45(twobody2,[0 25*day], [X0 0 0 V0], opts);

dist = sqrt(u45(:,1).^2+u45(:,2).^2);
vel = sqrt(u45(:,3).^2+u45(:,4).^2);
IT = find(dist<RT+RL,1);
fprintf('3) Impatto dopo %f giorni, velocita %f km/sec \n',t45(IT)/day,vel(IT));

figure(3); 
subplot(121); plot(u45(:,1),u45(:,2),'b');
hold on; patch(RT*sin(ang),RT*cos(ang),'r'); xlabel('x'); ylabel('y');
subplot(122); plot(t45,dist,'b',t45,(RT+RL)*ones(size(t45)),'r--'); xlabel('t'); ylabel('|r|');

