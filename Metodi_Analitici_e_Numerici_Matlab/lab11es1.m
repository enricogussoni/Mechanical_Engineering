clc; clear all; close all
%Es1 punto 1 (teorico)

%Es2 punto 2
x0 = 0;
xL = 1;
h = 0.1;


xnod = [x0:h:xL]';
xdof = xnod(2:end);
N = length(xdof);

% matrice
a = 1;    % diffusione
b = 1;    % trasporto
c = 1;    % reazione

U  = ones(N,1);
U1 = ones(N-1,1);

A = a/h  * ( 2*diag(U) - diag(U1,1) - diag(U1,-1) )+ ...
    b    * ( 1/2*diag(U1,1) - 1/2*diag(U1,-1) ) + ...
    c*h  * ( 2/3*diag(U) +1/6*diag(U1,1) + 1/6*diag(U1,-1) );

A(N,N) = a/h + b/2 + c*h/3;

frhs = @(x) 5*x.^3+15*x.^2-33*x-3;
F = h*frhs(xdof);
gN = 12;
F(N) = (h/2)*frhs(xdof(N)) + gN;

EIG = min(eig(A))    % la matrice e' definita positiva
if (A == A')
    disp ('matrice simmetrica')
else
    disp('matrice non simmetrica')
end% la matrice non e' simmetrica   

u = A\F;
u = [0;u]

plot(xnod,u,'bo-')

uex = @(x) 5*x.^3 - 3*x;
xdis = linspace(x0,xL,1000);
hold on; plot(xdis,uex(xdis),'k')
