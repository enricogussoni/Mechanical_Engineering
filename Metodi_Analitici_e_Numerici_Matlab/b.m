function [x, R, niter] = newtonsys(F, JF, x0, tol, nmax)

x=zeros(nmax);
R=1;
niter=1;
x(1)=x(0);

while (R>tol && niter<nmax)
    if JF==0...
        
       x(niter+1)=x(niter)-F(x(niter))\JF(x(niter));
       R=abs(x(niter+1)-x(niter));
       niter=niter+1;
end