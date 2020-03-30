% Elementi Finiti per problemi parabolici

% ut - uxx + ux + u = 5*x.^3 + 15*x.^2 - 33*x - 3 , x€[0,1] , t>0
% u(0,t)=0 t>0 % Condizione di Dirichlet a sinistra
% ux(1,t)=12 t>0 % Condizione di Neumann a destra
% u(x,0)=0 x€[0,1]

%% 1
h=0.2;
x0=0;
xL=1;
a=1; b=1; c=1;

xgdl=[x0+h:h:xL]';
n=length(xgdl);

dl=ones(1,n);
dc=ones(1,n-1);

Aa = 2*diag(dl) - diag(dc,1) - diag(dc,-1);
Aa(n,n) = 1;
Ab = -0.5*diag(dc,-1) +0.5*diag(dc,1);
Ab(n,n) = 1/2;
Ac = 2/3 * diag(dl) + 1/6 * diag(dc,1) + 1/6 * diag(dc,-1);
Ac(n,n) = 1/3;
A = a/h * Aa + b * Ab + c*h * Ac;

M=h*Ac;

gN=12;
f = @(x) 5*x.^3 + 15*x.^2 - 33*x - 3;
F = h*f(xgdl);
F(n)= h/2 * f(xgdl(n)) + gN;

xdef=[x0;xgdl];

%% 2
theta=0.5; % Crank-Nicolson
tau=0.01;
nit=30;
u=zeros(n,1);
udef=[0;u];

matrice = 1/tau .* M + theta .* A;
termine_noto = 1/tau .* M - (1-theta) .* A;

[L,U,P]=lu(matrice);

figure(1)
plot(xdef,udef)
title('Soluzione al variare di t')
hold on

for c=1:nit
    b=termine_noto*u + F;
    y=fwsub(L,P*b);
    u=bksub(U,y);
    udef=[0;u];
    plot(xdef,udef);
    drawnow; pause(.2)
end

%% 3
theta=0; % Eulero esplicito (in avanti)
u=zeros(n,1);
udef=[0;u];

matrice = 1/tau .* M + theta .* A;
termine_noto = 1/tau .* M - (1-theta) .* A;

[L,U,P]=lu(matrice);

figure(2)
plot(xdef,udef)
title('Soluzione al variare di t')
hold on

for c=1:nit
    b=termine_noto*u + F;
    y=fwsub(L,P*b);
    u=bksub(U,y);
    udef=[0;u];
    plot(xdef,udef);
    drawnow; pause(.2)
end

%% 4
theta=0; % Eulero esplicito (in avanti)
tau=0.001;
u=zeros(n,1);
udef=[0;u];

matrice = 1/tau .* M + theta .* A;
termine_noto = 1/tau .* M - (1-theta) .* A;

[L,U,P]=lu(matrice);

figure(3)
plot(xdef,udef)
title('Soluzione al variare di t')
hold on

for c=1:nit
    b=termine_noto*u + F;
    y=fwsub(L,P*b);
    u=bksub(U,y);
    udef=[0;u];
    plot(xdef,udef);
    drawnow; pause(.2)
end