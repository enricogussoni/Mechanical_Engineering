function [T,uh,vh] = leapfrog(fun,dfun,tmax,u0,v0,h)
% Il metodo richiede in input la funzione fun che descrive il problema di 
% Cauchy di ordine 2 definita come anonymous function (fun = @(t,u,v) ...), 
% la funzione dfun (anch'essa definita come anonymous function dfun = @(t,u,v)
% che contiene l'espressione di df/dy0, l'istante finale tmax dell'intervallo
% temporale di soluzione (l'istante iniziale è sempre t0 = 0), i dati iniziali 
% del problema di Cauchy u0 e v0 ed il passo di discretizzazione temporale h.
% Il metodo restituisce in output il vettore T degli istanti temporali e i 
% vettori u h e v h contenente gli spostamenti e le velocità calcolate numericamente.

T=0:h:tmax;
N=length(T);
uh=zeros(1,N);
vh=zeros(1,N);
uh(1)=u0;
vh(1)=v0;

tol=1e-12;
nmax=100;

for t=2:N
    
    uh(t)= uh(t-1) + h*vh(t-1) + ((h^2)/2)*fun(T(t-1),uh(t-1),vh(t-1));
    
    uv=uh(t-1);
    vv=vh(t-1);
    tv=T(t-1);
    un=uh(t);
    tn=T(t);
    
    % Funzioni per il metodo di Newton:
    f  = @(v) vv + (h/2) * (fun(tv,uv,vv) + fun(tn,un,v)) - v;
    df = @(v) (h/2)*dfun(tn,un,v) - 1;
    
    % Sottoiterazioni del metodo di Newton
    err = 1;
    k = 0;
    xv = uv;
    while (k< nmax && err> tol)
        dfx = df(xv);
        if dfx == 0
            error('Arresto per azzeramento di dfun');
        else
            xn = xv - f(xv)/dfx;
            err = abs(xn-xv);
            k = k+1;
            xv = xn;
        end
    end
    
    vh(t)=xn;
    
end    
