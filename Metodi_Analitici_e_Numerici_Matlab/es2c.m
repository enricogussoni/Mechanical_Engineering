x0 = [1;1;1]; 
tol = 1e-10; 
nmax = 10;

[x , R , niter ] = newtonsys(@F, @JF, x0, tol, nmax)
