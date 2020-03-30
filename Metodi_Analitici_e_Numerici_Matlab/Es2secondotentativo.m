%Es 2.2

%1
F = @(x) [x(1)^2 + x(2)^2 - 1; sin(x(1)*pi/2) + x(2)^3];

%2
JF= @(x) [2*x(1), 2*x(2); 0.5*pi*cos(x(1)*pi/2), 3*x(2)^2];

%4
x0=[1;1];
tol=1e-5;
nmax=10;

[F_si_annulla_per, residuo,num_iterazioni]=newtonsys(F,JF,x0,tol,nmax);