function [t_h,u_h,iter_nwt] = eulero_indietroGUS(f,df,t_max,y0,h)

u_h=zeros(t_max,1);
t_h=y0:h:t_max;
u(1)=y0;
N=length(t_h);

% parametri per mN
nmax = 100;
toll=1e-12;
iter_nwt=zeros(1,N);

for i=2:N
    u_old=u_h(it-1);
    t_new=t_h(it);
    % fun=df(x(i-1))*(x(i-1)-x(i));

% Da qui alla fine del ciclo while è inutile nel caso di funzioni lineari.
% In caso di f lineare in y si può ricavare direttamente.
    % funzioni per mN
    phi=@(u) u_old + h*f(t_new,u) - u;
    dphi=@(u) h*df(t_new,u) - 1;
    
    % sottoiterazioni per mN
    % (va implementata ogni volta da capo!)
    err=1;
    k=0;
    xv=u_old;
    while (k<nmax && err>toll)
    ... script di newton...
    end
%    u_h(it)=xn;
%    iter_nwt(it)=k;
end
