 %% ESERCIZIO 2

clc; clear all; close all


%Es2 punto 2
x0 = 0;
xL = 1;
h = 0.1;


xnod = [x0:h:xL]';
xdof = [x0:h:xL-h]';
N = length(xdof);

% matrice
a = 1;    % diffusione
c = 2;    % reazione

U  = ones(N,1);
U1 = ones(N-1,1);

A = a/h  * ( 2*diag(U) - diag(U1,1) - diag(U1,-1) )+ ...
    c*h  * ( 2/3*diag(U) +1/6*diag(U1,1) + 1/6*diag(U1,-1) );

A(1,1) =  a/h + c*h/3 + 1;

frhs = @(x) 2*x.^2-4*x;
F = h*frhs(xdof);
r = 3;
F(1) = (h/2)*frhs(xdof(1))+r;   %sovrascrivo la posizione 1 del vettore
 
%Es2 punto 3
if min(eig(A))>0
    fprintf('definita positiva \n')
else 
    fprintf('non definita positiva \n')
end

%Es 

if (A == A')
    fprintf('simmetrica \n')
else 
    fprintf('non simmetrica \n')
end
    

u = pcg(A,F)
u = [u; 0];

plot(xnod,u,'bo-')

uex = @(x) x.^2 -2*x +1;
xdis = linspace(x0,xL,1000);
hold on; plot(xdis,uex(xdis),'k')


