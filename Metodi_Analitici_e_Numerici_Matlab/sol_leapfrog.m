function [T,uh,vh]=sol_leapfrog(f,df,t_max,u0,v0,h)

% [T,u_h,v_h]=leapfrog(f,df,t_max,u0,v0,h)
% Risolve l'ode del secondo ordine:
%
% y'' = fun(t,y,y')
% y(0) = u0
% y'(0) = v0
%
% con il metodo leapfrog.
% Input:
% --> fun: e' la funzione che definisce il problema di Cauchy, dichiarata
%          come anonymous function:
%          fun = @(t,u,v) ...
% --> df : e' la funzione che definisce la derivata df/dy', dichiarate
%          come anonymous function:
%          df = @(t,u,v) ...
% --> t_max: l'istante finale dell' intervallo temporale di soluzione
%            (l'istante iniziale e' t_0=0)
% --> u0: il dato iniziale y(0) del problema di Cauchy
% --> v0: il dato iniziale y'(0) del problema di Cauchy
% --> h: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% --> T: vettore degli istanti in cui si calcola la soluzione discreta
% --> u_h: la soluzione discreta y calcolata negli istanti temporali
% discreti
% --> v_h: la soluzione discreta y' calcolata negli istanti temporali
% discreti

T=0:h:t_max;
uh = u0;
vh = v0;

for n=1:length(T)-1
    uh(n+1) = uh(n) + h*vh(n) + ((h^2)/2)*f(T(n),uh(n),vh(n));
    phi = @(w) vh(n) + (h/2)*( f(T(n),uh(n),vh(n)) + f(T(n+1),uh(n+1),w) ) - w;
    dphi = @(w) (h/2)*( df(T(n+1),uh(n+1),w) ) - 1;
    
    % Newton subiterations
    w = vh(n);
    wold = w;
    incr = 1; iter = 0;
    while incr > 1e-12 && iter < 100
        w = wold - phi(wold)/dphi(wold);
        incr = abs(w-wold);
        iter = iter + 1;
        wold = w;
    end
    vh(n+1) = w;
end
