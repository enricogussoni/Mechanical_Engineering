% Lab 12 Es1
% Elementi Finiti per problemi parabolici
% Problema di Cauchy-Dirichlet omogeneo per l'equazione di diffusione-trasporto-reazione

% ut - uxx + ux + u = 5*x.^3 + 15*x.^2 - 33*x - 3 , x€[0,1], t>0
% u(0,t)=0 , t>0
% ux(1,t)=12 , t>0
% u(x,0)=0 , x€[0,1]

%% 1
h=0.2;
x0=0;
xL=1;

xgdl = [x0+h:h:xL]';
n=length(xgdl);

a=1;
b=1;
c=1;

dl=ones(n,1);
dc=ones(n-1,1);

% matrice di massa M
M = h * (2/3 * diag(dl) + 1/6 * diag(dc,1) + 1/6 * diag(uc,-1));
M(n,n) = h * 1/3; % condizione di Neumann => 1/2 * h/3 (mezza campana)

% matrice di rigidezz A
Aa = a/h * (2 * diag(dl) - diag(dc,1) - diag(uc,-1));
Aa(n,n)= a/h * 1;

Ab = ...;
Ab(n,n)=b * 1/2;

Ac = ...;
Ac(n,n)=c*M;

A=Aa+Ab+Ac;

% termine noto
gN=12; % guarda i segni dalla derivazione per parti

f = @(x) 5*x.^3 + 15*x.^2 - 33*x - 3;
F = h*f(xgdl);
F(end) = h/2 * f(xgdl(end));

%% 2
dt=0.01;
npt=30;

u=zeros(n,1); % per salvare la soluzione e dato iniziale
xdef = [x0;xgdl];

figure(1)
plot(xdef,[0,u]);
hold on

theta = 0.5; % definisce il metodo di Crank-Nicholson

[L,U,P]=lu(M/dt + theta*A); 

for c=1:npt % ad ogni passo temporale cambia solo il termine noto
    b=(M/dt - (1-theta)*A)*u+F; % termine noto
    
    y=fwsub(L,P*b);
    u=bksub(U,y);
    
    plot(xbc,[0,u]);
    drawnow; 
    pause(.1) % per non plottare troppo veloce, solo comodità visiva
end
%  non sono richieste analisi dei dati, quindi è inutile memorizzare tutti
%  i risultati istante per istante: si possono sovrascrivere ad ogni
%  iterazione

%% 3
... copia-incolla con theta=0

%% 4
% bisogna controllare che il passo dt sia minore di 2 fratto il maggiore degli
% autovalori
EIG = eig(-A,M); % EIG=eig(-inv(M)*A)
2/max(EIG)

    