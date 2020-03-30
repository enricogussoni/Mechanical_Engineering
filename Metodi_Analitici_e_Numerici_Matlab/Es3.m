%Lab 8 es 3

GM=398600; % Km/s^2

%1 Orbita della Luna
% dw1 = w3;
% dw2 = w4;
% dw3 = -(GM/(w1^2 + w2^2)^(3/2))*w1;
% dw4 = -(GM/(w1^2 + w2^2)^(3/2))*w2;

%2
% x(0) = 405500; %Km
% y(0) = 0;
% dx(0) = 0;
% dy(0) = 0.964; %Km/s

sec_anno = 60*60*24*365 + 60*60*4;
duecorpi = @(t,y) [y(3); y(4); -GM*y(1)/(y(1)^2 + y(2)^2)^(3/2);...
                               -GM*y(2)/(y(1)^2 + y(2)^2)^(3/2)];
x0=405500;
v0=0.964;
[t45,u45]=ode45(duecorpi,[0 sec_anno],[x0 0 0 v0]);

num_timestep1=length(t45)

RT=6360; %Km
RL=1737; %Km
ang=[0:pi/100:2*pi];

figure(1)
subplot(121),
plot(u45(:,1),u45(:,2));
xlabel('x');
ylabel('y');
hold on;
patch(RT*sin(ang),RT*cos(ang),'r');
subplot(122);
plot(t45,sqrt(u45(:,1).^2+u45(:,2).^2)/x0);
xlabel('t');
ylabel('|r|');

%3
options=odeset('reltol',1e-6); % modifica tolleranza
[T45,Y45]=ode45(duecorpi,[0 sec_anno],[x0 0 0 v0],options);

num_timestep2=length(T45)

figure(2)
subplot(121),
plot(Y45(:,1),Y45(:,2));
xlabel('x');
ylabel('y');
hold on;
patch(RT*sin(ang),RT*cos(ang),'r');
subplot(122);
plot(T45,sqrt(Y45(:,1).^2+Y45(:,2).^2)/x0);
xlabel('t');
ylabel('|r|');

%4 Con attrito
C=1e-6;
duecorpi_attr = @(t,y) [y(3); y(4); -GM*y(1)/(y(1)^2 + y(2)^2)^(3/2) - C*y(3);...
                                    -GM*y(2)/(y(1)^2 + y(2)^2)^(3/2) - C*y(4)];

options=odeset('reltol',1e-6);
[T45a,Y45a]=ode45(duecorpi_attr,[0 25*24*60*60],[x0 0 0 v0],options);

num_timestep3=length(T45a)

figure(3)
subplot(121),
plot(Y45a(:,1),Y45a(:,2));
xlabel('x');
ylabel('y');
hold on;
patch(RT*sin(ang),RT*cos(ang),'r');
subplot(122);
plot(T45a,dist,T45a,(RT+RL)*ones(size(T45a)));
xlabel('t');
ylabel('|r|');

dist=sqrt(Y45a(:,1).^2 + Y45a(:,2).^2);
vel=sqrt(Y45a(:,3).^2 + Y45a(:,4).^2);
II=find(dist<RT+RL,1); % Istante di impatto (nel vettore dei tempi)
fprintf('Impatto dopo %f giorni a velocità %f km/sec \n',T45a(II)/(60*60*24),vel(II));