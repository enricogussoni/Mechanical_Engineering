F = @(x) [ x(1)^2 + x(2)^2-1; sin(pi*x(1)/2) + x(2)^3];

JF = @(x) [2*x(1) , 2*x(2); 0.5*pi*cos(0.5*pi*x(1)), 3*x(2)^2];
       
x0 = [1;1]; 
tol = 1e-5; 
nmax = 10;

[x , R , niter ] = newtonsys(F, JF, x0, tol, nmax);
