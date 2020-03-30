clc; clear all; close all
x0 = 0;
xL = 1;
h = 0.2;
xdof = [x0+h:h:xL]';
N = length(xdof);

a  = 1; b = 1; c = 1;
U  = ones(N,1);
U1 = ones(N-1,1);
% matrice di massa M
M = h *  ( 2/3*diag(U) +1/6*diag(U1,1) + 1/6*diag(U1,-1) );
M(N,N)= h*(1/3);
% matrice di rigidezza A
Aa = a/h * ( 2*diag(U) - diag(U1,1) - diag(U1,-1) );
Aa(N,N)=(a/h)*1;
Ab = b * ( 1/2*diag(U1,1) - 1/2*diag(U1,-1) );
Ab(N,N)=b*0.5;
Ac = c*M;
A  = Aa+Ab+Ac;
% termine noto
gN = 12;
frhs = @(x) 5*x.^3+15*x.^2-33*x-3;
F = h*frhs(xdof);
F(N) = (h/2)*frhs(xdof(N)) + gN;

u = zeros(N,1);
dt = 0.01;
xbc = [x0;xdof];
plot(xbc,[0;u],'ro-'); hold on

theta = 0.5;
[L,U,P]=lu(M/dt+theta*A);
for i=1:30
    b=(M/dt-(1-theta)*A)*u+F;
    y=fwsub(L,P*b);
    u=bksub(U,y);
    plot(xbc,[0;u],'bo-');
    drawnow; pause(.1)
end

EIG = eig(-A,M)
EIG2=eig(-inv(M)*A)
2/(max(abs(EIG)))