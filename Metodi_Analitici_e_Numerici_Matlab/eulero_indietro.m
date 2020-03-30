function [t_h,u_h,iter_nwt]=eulero_indietro(f,df,t_max,y_0,h)

% [t_h,u_h,iter_pf]=eulero_indietro(f,df,t_max,y_0,h)
% Risolve il problema di Cauchy
%
% y'=f(y,t)
% y(0)=y_0
%
% utilizzando il metodo di Eulero Indietro.
% Per ciascun istante temporale si calcola u_(n+1)
% trovando lo zero dell'equazione (a priori non lineare)
% F(u) = u_n + h*f(t_{n+1},u) - u = 0
% con il metodo di Newton. Per questo e' necessario passare
% in input l'espressione di df/dy, in modo da creare
% dF(u)/du = h* df/dy(t_{n+1},u) - 1.
%
% Input:
% -> f: function che descrive il problema di Cauchy
%       (dichiarata come anonymous function)
%       deve ricevere in ingresso due argomenti: f = @(t,y)
% -> df: function che descrive df/dy, anch'essa anonymous function
%        che riceve in ingresso 2 argomenti: df = @(t,y)
% -> t_max: l'istante finale dell' intervallo temporale di soluzione
%                 (l'istante iniziale e' t_0=0)
% -> y_0: il dato iniziale del problema di Cauchy
% -> h: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% -> t_h: vettore degli istanti in cui si calcola la soluzione discreta
% -> u_h: la soluzione discreta calcolata nei nodi temporali t_h
% -> iter_nwt: vettore che contiene il numero di iterazioni
%                del metodo di Newton necessarie a risolvere l'equazione
%                non lineare ad ogni istante temporale.

t0 = 0;
t_h = t0:h:t_max;
% inizializzo il vettore che conterra' la soluzione discreta
N   = length(t_h);
u_h = zeros(1,N);
u_h(1) = y_0;
% parametri per il metodo di Newton
nmax = 100;
toll = 1e-12;
iter_nwt=zeros(1,N);  %vettore in cui salvo per ogni passo temporale il numero di iterazioni effettuate dal metodo di Newton 

for it=2:N
    u_old = u_h(it-1);
    t_new = t_h(it);
    % Funzioni per il metodo di Newton:
    phi  = @(u) u_old + h * f( t_new, u ) - u;
    dphi = @(u) h*df(t_new, u) - 1;
    % Sottoiterazioni del metodo di Newton
    err = 1;
    k = 0;
    xv = u_old;
    while (k< nmax && err> toll)
        dfx = dphi(xv);
        if dfx == 0
            error('Arresto per azzeramento di dfun');
        else
            xn = xv - phi(xv)/dfx;
            err = abs(xn-xv);
            k = k+1;
            xv = xn;
        end
    end
    u_h(it) = xn;
    iter_nwt(it) = k;
end