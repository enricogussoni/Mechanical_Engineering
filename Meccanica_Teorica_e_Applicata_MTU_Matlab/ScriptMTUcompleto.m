clear all
close all
clc

%Carico dati (curva potenza/giri al minuto)
%Salvare script nella stessa cartella in cui ci sono le funzioni)
%Controllare di aver selezionato il path giusto nell'interfaccia di
%comando)

load Wmax.mat
load Wmin.mat

%Calcolo Coppia Massima e Minima

T_max_endo=W_max_endo./(rpm_endo*2*pi/60);
T_min_endo=W_min_endo./(rpm_endo*2*pi/60);

figure
plot(rpm_endo,W_max_endo/1000)
grid
xlabel('[rpm]')
ylabel('[kW]')

figure
plot(rpm_endo,T_max_endo)
grid
xlabel('[rpm]')
ylabel('[Nm]')
lim=axis;
xlim([0 lim(2)])
ylim([0 lim(4)])

%Costruzione fascio di rette per trovare il regime di coppia massima

m=[2:10:102];
y=2*pi/60*[0;rpm_endo]*m;

figure
plot(rpm_endo,W_max_endo/1000)
grid
hold on
plot([0;rpm_endo],y/1000) %si aggiunge lo zero per avere l'origine

%Curva di coppia la variare del grado di ammettenza (gamma, "apertura gas")

gamma=[0:.2:1];
T_endo=zeros(length(T_max_endo),length(gamma));
legenda=[];
figure

for ii=1:length(gamma)
    T_endo(:,ii)=gamma(ii)*T_max_endo+(1-gamma(ii))*T_min_endo;
    legend=[legenda;['\gamma=' sprintf('%4.2f',gamma(ii))]]; 
    %'\' serve per le lettere greche
    %per sprintf vedi help sprintf


hold on
plot(rpm_endo,T_endo)
grid
xlabel('[Nm]');
ylabel('[]');

end


%Dati autovettura

m=(1025+160+50+200); % kg
Cx=0.38;
S=2.5; % m^2
rho=1.2258; % kg/m^3
fv=0.01;
g=9.81; % m/s^2
Rr=0.27; % m (raggio delle ruote)


%Curva di carico per veicolo in piano (alpha=0)

alpha=-0.05;
vel=[0:190]; %vettore di velocità possibili (passo unitario)
Cr=(m*g*(sin(alpha)+fv*cos(alpha))+0.5*Cx*S*rho*(vel/3.6).^2)*Rr;

figure
%hold on
plot(vel,Cr)
xlabel('[Km/h]')
ylabel('[Nm]')
%title('Curva di carico a regime')
grid


%Coppia alle ruote %slide 10->

taud=1/4.071; %differenziale
tauc=[1/3.909 1/2.158 1/1.480 1/1.121 1/.921]*taud; %cambio
etac=[0.9,0.94,0.94,0.97,0.97]; 

v_endo=zeros(length(rpm_endo),length(tauc));
%%per gamma=1 si usa T_max_endo; 
%si considera quindi la coppia massima, altrimenti si sceglie
%T_endo(:,colonna della gamma scelta)

for jj=1:length(tauc)
    T_endo_r_g1(:,jj)=etac(jj)/tauc(jj).*T_endo(:,3); %jjesima marcia
    v_endo(:,jj)=(2*pi/60*rpm_endo*tauc(jj))*Rr*3.6; %km/h    
end

%figure
hold on
plot(v_endo,T_endo_r_g1)
xlabel('[km/h]')
ylabel('[Nm]')
grid

%Forza di traino residua

%%gamma=1;
v_endo=zeros(length(rpm_endo),length(tauc));

for jj=1:length(tauc)
    v_endo(:,jj)=(2*pi/60*rpm_endo*tauc(jj))*Rr*3.6; %km/h
    Faer=0.5*rho*Cx*S*(v_endo(:,jj)/3.6).^2;
    Ftraino(:,jj)=(T_endo_r_g1(:,jj)/Rr-Faer-m*g*fv);
end

figure
plot(v_endo,Ftraino)
xlabel('[km/h]')
ylabel('[N]')
grid

