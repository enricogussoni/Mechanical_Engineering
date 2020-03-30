clear all
clc

% Calcoli per un ventilatore assiale da airboat

%% Dati generici

g = 9.81;                % accelerazione di gravità [m/s^2]
visc_cin_a = 1.50e-5;    % viscosità cinematica dell'aria [m^2/s]
p_atm = 1e5;             % pressione atmosferica [Pa]

v = 100/3.6;             % velocità barchino [m/s] (100km/h)
T = 8000;                % spinta [N]
D = 2;                   % diametro del rotore [m]
d = 0.4;                 % diametro ogiva [m]
Dm=(D+d)/2;              % diametro medio [m]
b = (D-d)/2;             % altezza di pala
A = (D^2-d^2)/4 * pi;    % area di passaggio
ro_a = 1.225;            % densità dell'aria [kg/m^3]
m = A*v*ro_a;            % portata massica [kg/s]
u_tip = 100;             % velocità di trascinamento (per ip. non comprimibile) [m/s]
omega = u_tip/(D/2);     % velocità angolare delle pale
u = omega*(Dm/2);        % velocità di trascinamento mid della pala [m/s]

%% Lavoro e triangoli di velocità

v1t = 0;                 % velocità assoluta tangenziale in ingresso
v1 = v;                  % " " in ingresso
v1a=v1;                  % velocità assoluta assiale in ingresso
alpha1=pi/2;
w1 = sqrt(v1^2+u^2);     % velocità relativa in ingresso (nell'ipotesi di v1t=0)

l = v*T/m;               % lavoro specifico [J/kg] (l_eu - lw)
v2t= l/u;                % " " tangenziale " "
v2a= v1a;                % " " assiale " "
v2= sqrt(v2t^2+v2a^2);   % velocita assoluta in uscita [m/s]

alpha2= acos(v2t/v2);
alpha2deg= acos(v2t/v2)*180/pi;
w2a= v2a;                % " relativa " " "
w2t= v2t - u;            % " " tangenziale " "
w2 = sqrt(w2a^2+w2t^2);  % " " in uscita
beta2=acos(w2t/w2);
beta2deg=acos(w2t/w2)*180/pi;
beta1=atan(v1/u);
beta1deg=atan(v1/u)*180/pi;

% Verifica Criterio di De Haller
% if (w2/w1>0.72)
%     disp('Il criterio di De Haller è verificato')
% else
%     disp('Il criterio di De Haller non è verificato')
% end

pt1=p_atm+v1^2/2*ro_a;
pt2=pt1+ro_a*l;
p2=pt2-ro_a*v2^2/2;
X=(p2-p_atm)/(pt2-pt1);

%% Parametri adimensionali
lw_i=0;                              % lw di primo tentativo
H=(l-lw_i)/g;
omegas = omega*(A*v)^0.5/(g*H)^(3/4);
Ds = D*(g*H)^(1/4)/v^0.5;

% % numero di pale (basandosi sulle tabelle in funzione dell'ws)
% if (omegas>0 && omegas<2)
%     z=16;
% elseif(omegas>2 && omegas<5)
%     z=12;
% elseif (omegas>5 && omegas<6)
%     z=6;
% else fprintf('ws non accettabile');
% end

z = 6;                                 % numero di pale ipotizzato(da analisi dei modelli esistenti e su iterazioni del codice)                
s = 2*pi/z*Dm;                         % passo [m]

Re = (u/2)*(D-d)/2/visc_cin_a;         % Reynolds al tip
u_Re = [0:0.1:(D/2 - d/2)].*omega;     % vettore delle velocità lungo la pala
vet_Re = (u_Re./2)*(D-d)/2/visc_cin_a; % vettore dei numeri di Re lungo la pala

%% Dimensionamento del motore

P = l*m;                    % potenza resistente [W]
Pcv = P*1.36e-3;            % " in [CV]
C = P/omega;                % coppia all'asse [Nm]

% verifica a fatica dell'albero motore

%% Scelta del profilo

c = 0.45;                                   % lunghezza della corda [m]
solidity=c/s
xa=0;                                       % ascissa di max curvatura iniziale
mc=0.18;                                    % coefficiente della formula della deflessione cinematica (funzione di beta1)
bl=0.93;                                    % altro coefficiente di deflessione cinematica
camber=abs((pi-beta1)-beta2)*180/pi+0.5;    % deflessione Geometrica [deg]
AR= 2;                                      % Aspect Ratio ipotizzato (Danile Pizzo)
lc= b/AR;                                   % lunghezza corda (verifica)

sf=0.08;                                    % altezza sezione frontale singola pala
S=sf*b;                                     % sezione frontale singola pala

% DF = 1-w2/w1+(v2t-v1t)/(2*w1*solidity);    % Fattore di diffusione globale

i10= 0.5;                                   % incidenza i0_10 [deg] (da grafico Lieblein) (funzione di beta1)
delta10= 0.2;                               % deflessione cinematica delta0_10 [deg] (funzione di beta1)
nl=-0.11;                                   % coefficiente n di Lieblein (da grafico) (funzione di beta1)
i = i10+nl*camber;                          % incidenza ottimale
% deflessione = 30;                         % [deg] (da grafico in funzione di i)
cl= -0.15;                                  % coefficiente di lift da grafico per NACA 65-206
cd= 0.007;                                  % coefficiente di drag " " " " "

delta=delta10+mc*camber*(solidity^(-bl));   % deflessione cinematica [deg]

gamma=(90-(beta1*180/pi+i));                % stagger (angolo di calettamento) [deg]

lw=z*0.5*cd*v^3*ro_a*S/m;                   % perdita sui profili [J/kg]
rendimento= (l-lw)/l;                       

% Definizione: c_p = (pt2-pt1)/(pt1-p1) con pt=pressione totale
% c_p = 0.03;                                 % coefficiente di perdita (da grafico in funzione di i)
 
%% Calcolo alternativo del numero di pale
% da tesi di Daniele Pizzo

% t=lc/(solidity*cos(gamma));               % passo palare (Howell)
% zc=(pi*(D+d)/2)/t;                        % numero di pale calcolato (Howell)
% n_pale=ceil(zc);                          % numero di pale effettivo (Howell)

%% Svergolamento dell'elica
% Progetto del rotore a vortice libero

Dm=(D+d)/2;                                  % diametro medio
v_ang=u/(Dm/2);                              % velocità angolare

% r = linspace(d/2,D/2,10);                    % vettore di quote radiali
% u_r = v_ang*r;                               % distribuzione di velocità di trascinamento
K = v2t*Dm/2;                                % costante Vt*r

% Controllo sulla non espansione alla base
% if (X(1)<0)
%     disp('Il grado di reazione alla base della pala è negativo')
% end

%% HUB

uH=d/2*omega;
v2tH= l/uH;                 % " " tangenziale " "
v2aH= v1a;                  % " " assiale " "
v2H= sqrt(v2tH^2+v2aH^2);   % velocita assoluta in uscita [m/s]

alpha2H= acos(v2tH/v2H);
alpha2degH= acos(v2tH/v2H)*180/pi;
w2aH= v2aH;                 % " relativa " " "
w2tH= v2tH - uH;            % " " tangenziale " "
w2H = sqrt(w2aH^2+w2tH^2);  % " " in uscita
beta2H=acos(w2tH/w2H);
beta2degH=acos(w2tH/w2H)*180/pi;
beta1H=atan(v1/uH);
beta1degH=atan(v1/uH)*180/pi;

camberH=abs((pi-beta1H)-beta2H)*180/pi+0.5;    % deflessione Geometrica [deg]

i10H= 2;                                       % incidenza i0_10 [deg] (da grafico Lieblein) (funzione di beta1)
delta10H= 1;                                   % deflessione cinematica delta0_10 [deg] (funzione di beta1)
nlH=-0.47;                                     % coefficiente n di Lieblein (da grafico) (funzione di beta1)
iH = i10H+nlH*camberH;                         % incidenza ottimale [deg]

%% TIP

uT=100;
v2tT= l/uT;                 % " " tangenziale " "
v2aT= v1a;                  % " " assiale " "
v2T= sqrt(v2tT^2+v2aT^2);   % velocita assoluta in uscita [m/s]

alpha2T= acos(v2tT/v2T);
alpha2degT= acos(v2tT/v2T)*180/pi;
w2aT= v2aT;                 % " relativa " " "
w2tT= v2tT - uT;            % " " tangenziale " "
w2T = sqrt(w2aT^2+w2tT^2);  % " " in uscita
beta2T=acos(w2tT/w2T);
beta2degT=acos(w2tT/w2T)*180/pi;
beta1T=atan(v1/uT);
beta1degT=atan(v1/uT)*180/pi;

camberT=abs((pi-beta1T)-beta2T)*180/pi+0.5;    % deflessione Geometrica [deg]

i10T= 0.5;                                     % incidenza i0_10 [deg] (da grafico Lieblein) (funzione di beta1)
delta10T= 0.12;                                % deflessione cinematica delta0_10 [deg] (funzione di beta1)
nlT=-0.11;                                     % coefficiente n di Lieblein (da grafico) (funzione di beta1)
iT = i10T+nlT*camberT;                         % incidenza ottimale
