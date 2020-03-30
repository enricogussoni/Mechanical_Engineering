clear all
close all
clc
warning off

global G
global a s b c 
global m1 J1 m2 J2 m3 J3 k1 r1 k2 r2
global z0 wz psiz
global theta0 g0 dg0 ddg0 delta0 ddelta0 dddelta0 DL01 theta10 theta20 x p
global g_old delta_old
global modi_i modi_i_smo

% dati
G = 9.81;           % gravità [m/s^2] 
a = 2.520;          % posizione baricentro G1 [m]
x = 2;       
p = 2;       
b = 7.493;          % posizione baricentro G2 [m]
c = 9.02;           % posizione baricentro G3 e punto d'applicazione della forzante
s = 3.20;          
m1 = 2280;          % massa braccio 1 [kg]
J1 = m1*(a+s)^2/12; % momento inerzia braccio 1 [kg*m^2]
m2 = 870;           % massa braccio 2 [kg]
J2 = m2*(p+x)^2/12; % momento inerzia braccio 2 [kg*m^2]
m3 = 1120;          % massa martello (concentrata) [Kg]
J3 = 0;             % momento inerzia martello [kg*m^2]
k1 = 2*18096000;    % rigidezza molla 1 [N/m]
k2 = k1/2;          % rigidezza molla 2 [N/m]
r1 = 2*50000;         % smorzatore 1 [Kg/s]
r2 = r1/2;          % smorzatore 2 [Kg/s]

% forzamento
z0 = 137000;        % ampiezza [N] (da D_martello=0.118 m e P_martello=125 bar)
wz = 80;            % pulsazione (da f=800 colpi/min)         
psiz = 0;           % fase

%% 1 Grado di Libertà

% equilibrio statico
theta0 = pi/4;      % angolo di equilibrio
g_old = 2.8;
delta_old = pi/3;

[g0,dg0,ddg0,delta0,ddelta0,dddelta0] = cinematica1(theta0);
DL01 = precarico1();

% dati di simulazione
x0 = [0; theta0]; % condizioni iniziali
T = 10;                 % [s]   durata
dt = 5e-2;              % [s]   passo

% integrazione del sistema non lineare
[t_nl,y_nl] = rk4('manovellismo_non_lineare',T,dt,x0);

% integrazione del sistema linearizzato
[t_l,y_l] = rk4('manovellismo_lineare',T,dt,x0);

% confronto tra sistema linearizzato e non linearizzato
figure(1)
plot(t_nl,y_nl(2,:)*180/pi,t_l,y_l(2,:)*180/pi,'g');
ylabel('Rotazione del primo braccio [deg]')
xlabel('Tempo [s]')
title ('risposta nel tempo forzata 1 gdl')
legend('Non lineare','Linearizzato')
%%

ddy1_0 = -a*sin(theta0);
ddy2_0 = -(a+s)*sin(theta0)-p*sin(theta0-(pi*10)/24);
ddy3_0 = -(a+s)*sin(theta0)-(p+x)*sin(theta0-(pi*10)/24);

M0 = energia_cinetica();
K0 = k1*(DL01*ddg0+dg0^2)+G*(m1*ddy1_0+m2*ddy2_0+m3*ddy3_0); 

w0 = sqrt(K0/M0)
f0 = w0/(2*pi);

R0 = r1*dg0^2;                                              
rc = 2*M0*w0;
h = R0/rc
wd = w0*sqrt(1-h^2);
fd = wd/(2*pi);



%% 2 Gradi di Libertà

theta10 = theta0;
theta20 = pi/3;

b2 = sqrt(((a+s)*cos(theta0)+p*cos(theta0-(pi*10)/24))^2+((a+s)*sin(theta0)+p*sin(theta0-(pi*10)/24))^2);
b3 = sqrt(((a+s)*cos(theta0)+(p+x)*cos(theta0-(pi*10)/24))^2+((a+s)*sin(theta0)+(p+x)*sin(theta0-(pi*10)/24))^2);

MAS = [m1,m2,m3];
mfis = diag(MAS);
KAP = [k1,k2];
Kfis = diag(KAP);
SMORZ = [r1,r2];
Rfis = diag(SMORZ);

DL02 = precarico2();

[g,dg,ddg,delta,ddelta,dddelta] = cinematica1(theta0);
[k,dk,ddk,tau,dtau,ddtau] = cinematica2(theta20);

M = [m1*a^2+J1+m2*b2^2+J2+m3*b3^2+J3,m2*p*b2+J2+m3*b3*(p+x)+J3;m2*p*b2+J2+m3*b3*(p+x)+J3,m2*p^2+J2+(p+x)^2*m3+J3];
K1 = [dg 0;0 dk]'*Kfis*[dg 0;0 dk];
R = [dg 0;0 dk]'*Rfis*[dg 0;0 dk];  
K2 = [k1*ddg*DL01 0;0 k2*ddk*DL02];
H1 = [-a*sin(theta10) 0;0 0];
H2 = [-(a+s)*sin(theta10)-p*sin(theta10+theta20+7/6*pi) -p*sin(theta10+theta20+7/6*pi);-p*sin(theta10+theta20+7/6*pi) -p*sin(theta10+theta20+7/6*pi)];
H3 = [-(a+s)*sin(theta10)-(p+x)*sin(theta10+theta20+7/6*pi) -(p+x)*sin(theta10+theta20+7/6*pi);-(p+x)*sin(theta10+theta20+7/6*pi) -(p+x)*sin(theta10+theta20+7/6*pi)];
K3 = G*(m1*H1+m2*H2+m3*H3);
K = K1+K2+K3;
B = z0*[c*cos(theta10-pi/4),(p+x)*cos(theta10+theta20+7/6*pi)]';          

% No smorzamento
[modi_i,auto_val] = eig(M\K);   
omega_i = sqrt(diag(auto_val))          % frequenze proprie
modi_i                                  % modi di vibrare
auto_val

% Si smorzamento
A = [M,zeros(2,2);zeros(2,2),M]\[-R,-K;M,zeros(2,2)]; 
[modi_i_smo,auto_val_smo] = eig(A);  
omega_i_smo = sqrt(diag(auto_val_smo))  % frequenze proprie
modi_i_smo                              % modi di vibrare
auto_val_smo

% No smorzamento
freq = (0:.01:10);
for ii = 1:length(freq)
Omega(ii) = 2*pi*freq(ii);
C1(ii,:) = inv(-(Omega(ii)^2)*M + K)*B;    
end

figure(2)
plot(Omega,abs(C1));
ylabel('Ampiezza relativa oscillazione x0/F0???')
xlabel('Frequenze [Hz]')
title ('Risposta in frequenza 2 gdl non smorzata')
legend('theta 1','theta 2')
figure(3)
plot(Omega,angle(C1));
ylabel('fase oscillazione')
xlabel('Frequenze [Hz]')
title ('Risposta in frequenza 2 gdl non smorzata')
legend('theta 1','theta 2')


% Si smorzamento
for ii = 1:length(freq)
Omega(ii) = 2*pi*freq(ii);
C2(ii,:) = inv(-(Omega(ii)^2)*M + 1i*Omega(ii)*R + K)*B;    
end

figure(4)
plot(Omega,abs(C2));
ylabel('Ampiezza relativa oscillazione x0/F0???')
xlabel('Frequenze [Hz]')
title ('Risposta in frequenza 2 gdl smorzata')
legend('theta 1','theta 2')
hold off;
figure(5)
plot(Omega,angle(C2));
ylabel('fase oscillazione')
xlabel('Frequenze [Hz]')
title ('Risposta in frequenza 2 gdl smorzata')
legend('theta 1','theta 2')

% Risposta nel tempo

Q= [0,0]               %risposta nel tempo non forzata
dz= @(t,y) [y(3);y(4);...
    (Q(1)-((M(1,2)/M(2,2))*(Q(2)-R(2,1)*y(3)-R(2,2)*y(4)-K(2,1)*y(1)-K(2,2)*y(2))+R(1,1)*y(3)+R(1,2)*y(4)+K(1,1)*y(1)+K(1,2)*y(2)))/(M(1,1)-(M(1,2)*M(2,1)/M(2,2)));...
    (Q(2)-(M(2,1)*((Q(1)-((M(1,2)/M(2,2))*(Q(2)-R(2,1)*y(3)-R(2,2)*y(4)-K(2,1)*y(1)-K(2,2)*y(2))+R(1,1)*y(3)+R(1,2)*y(4)+K(1,1)*y(1)+K(1,2)*y(2)))/(M(1,1)-(M(1,2)*M(2,1)/M(2,2))))+R(2,1)*y(3)+R(2,2)*y(4)+K(2,1)*y(1)+K(2,2)*y(2)))/M(2,2)]
   inizio = [0.01;0.01;0;0];      % perturbazione iniziale tattico fallica oscar e w 5894
[tode, yode]=ode45(dz,0:5e-3:100,inizio);


figure

plot (tode, (theta10 + yode(:,1))*180/pi,'b')
hold on
grid on
plot (tode, (theta20 + yode(:,2))*180/pi,'r')
ylabel('Rotazione dei bracci [deg]')
xlabel('Tempo [s]')
title ('risposta nel tempo non forzata 2 gdl')
legend('theta 1','theta 2')
hold off


% 1- diminuire omega per visualizzare il contributo del forzamento nell'1
% gdl 2-riportare lo smorzamento al valore corretto 4- RISPOSTA FORZATA NEL TEMPO 2 GDL!!!!!
% (Gus: nella consegna del progetto che ha caricato Egidio su beep non c'è
% scritto di farla, eppur va fatta)
% 5- RISPOSTA NON FORZATA 1 GDL % 6- Controllare i Labels dei grafici 

% Dal file di Egidio:
% PROGETTO D’ANNO - MECCANICA DELLE VIBRAZIONI

% Lo svolgimento del progetto deve essere realizzato in gruppi di massimo 4 componenti ciascuno. 
% Tale svolgimento sarà descritto in un elaborato da consegnarsi prima dell’inizio della sessione d’esame in data che verrà comunicata dal docente. La valutazione del progetto sarà la medesima per tutti i componenti del gruppo.
% In sede di esame orale verrà discusso il progetto e potranno essere visionati i codici Matlab utilizzati per lo svolgimento. La discussione del progetto sarà parte integrante dell’esame orale e sarà pertanto individuale.
% Il progetto si articolerà su tre diversi punti:
% 
% Punto 1: Individuare un sistema vibrante reale e impostarne una schematizzazione volta ad identificare:
% Masse in movimento con rappresentazione dei vincoli e dei legami cinematici FATTO
% Elementi elastici e smorzanti FATTO
% Forzanti FATTO
% 
% Punto 2: Semplificare il sistema fino a schematizzarlo come un sistema vibrante non lineare ad 1 gdl e risolverlo, mediante l’ausilio di codici Matlab. 
% Si calcolino:
% la risposta del sistema non lineare: sistema forzato e moto libero assegnate condizioni iniziali; FATTO
% la risposta del sistema linearizzato nell’intorno di una posizione di equilibrio: sistema forzato e moto libero assegnate condizioni iniziali; FATTO
% la frequenza propria del sistema linearizzato. FATTO

% I risultati ottenuti dovranno essere commentati e, eventualmente, valutati al variare di alcuni variabili.
% 
% Punto 3: Semplificare il sistema fino a schematizzarlo come un sistema vibrante linearizzato a 2 gdl e risolverlo mediante l’ausilio di codici Matlab. Si calcolino:
% le frequenze proprie e i modi di vibrare del sistema in assenza e in presenza di smorzamento; FATTO
% le risposte del sistema libero in assenza e in presenza di smorzamento assegnate le condizioni iniziali; FATTO
% le funzioni di risposta in frequenza tra l’ingresso e le uscite del sistema. FATTO

% I risultati ottenuti dovranno essere commentati e valutati, eventualmente, per diversi livelli di smorzamento.