clear all
clc
close all

warning off %disabilita i warning che non interrompono l'esecuzione
%% Dati
m=2000; %kg
k=1000; %N/m
x0=0; %m
xp0=10; %m/s

%% Moto Sistema Vibrante smorzato con h<1

h=0.05;
ome0=sqrt(k/m);
alpha=h*ome0;
omed=ome0*sqrt(1-h^2);

A=x0;
B=(xp0+alpha*x0)/omed;

t=[0:0.01:30];

x=exp(-alpha*t).*(A*cos(omed*t)+B*sin(omed*t));

figure
plot(t,x)
grid
ylabel('spostamento [m]');
xlabel('tempo [s]');

%% Moto Sistema Vibrante smorzato con h=1

h=1;
alpha=h*ome0;

X1=x0;
X2=xp0+alpha*x0;

x2=exp(-alpha*t)*X1+X2*t.*exp(-alpha*t);

hold on %per plottare sullo stesso grafico di prima
plot(t,x2,'r')

% stesse condizioni iniziali, massimo spostato a sinistra, ampiezza minore,
% minor tempo di riequilibrio

%% %% Moto Sistema Vibrante smorzato con h>1

h=1.5;
alpha=h*ome0;

alpha1=-alpha+sqrt(h^2-1)*ome0;
alpha2=-alpha-sqrt(h^2-1)*ome0;

X1=(-xp0+alpha2*x0)/(alpha2-alpha1);
X2=(xp0-alpha1*x0)/(alpha2-alpha1); 

% anche se uso ancora i nomi X1 e X2 per le variabili questo non cambia il
% grafico gi� fatto (sequenzialit� di Matlab), al pi� si perdono i valori
% precedenti in caso di salvataggio

X3=X1*exp(alpha1*t)+X2*exp(alpha2*t);

hold on
plot(t,X3,'g')


%% Moto forzato (a regime)
m=20; %kg
k=19; %N/m

a=[0:0.05:3];
h=[0.05 0.1 0.2 0.3 1/sqrt(2) 1.5];

figure
hold all %mantiene tutti i grafici

for kk=1:length(h)
    XsuXst=1./sqrt((1-a.^2).^2+(2*a.*h(kk)).^2); %modulo
    plot(a,XsuXst)
end

grid
ylabel('a')
xlabel('|X/Xst|')

figure %evita di coprire il grafico precedente
hold hall %mantiene tutti i grafici

for kk=1:lenght(h)
    RXsuXst=atan2d(-2*a*h(kk),1-a.^2); %fase (il comando atan2d mette numeratore,divisore 
                                       %invece di richiedere / e parentesi)
    plot(a,XsuXst)
end

grid
ylabel('a')
xlabel('R|X/Xst|')

%% Moto forzato (con transitorio -> risposta completa)

m=2000;
r=2513;
k=29740;

ome0=sqrt(k/m);
omed=ome0*sqrt(1-h^2);
h=r/(2*m*ome0);
alpha=r/(2*m);

t=[0:0.01:30];

F0=10000; %N
omef=10; %[rad]/s, pulsazione della forzante 
psi=0; %[rad], sfasamento forzante

% condizioni iniziali
x0=0;  
xp0=10;

frf_F=1/(-omef^2*m+j*omef*r+k);
mod_F=abs(frf_F);
phi_F=angle(frf_F);

F=F0*cos(omef*t+psi);

xsig=F0*mod_F;

A1=x0-xsig*cos(psi+phi_F); %lo spostamento � sfasato rispetto alla forzante
                          %-> sfasamento relativo
B2=(alpha*x0+xp0-alpha*xsig*cos(psi+phi_F)+xsig*omef*sin(psi+phi_F))/omed;

x=exp(-alpha*t).*(A1+cos(omed+t)+B2*sin(omed*t))+xsig*cos(omef*t+psi+phi_F);

figure
plotyy(t,F/1000,t,x) %plotta sullo stesso diagramma grandezze con scale diverse

% A regime sono in zona quasistatica (forzante e spostamento a regime)con
% omef=1 e in zona sismografica con omef=10