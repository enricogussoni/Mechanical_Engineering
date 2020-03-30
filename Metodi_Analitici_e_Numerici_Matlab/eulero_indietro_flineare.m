function [t_h,u_h]=eulero_indietro_flineare(t_max,y_0,h)

% [t_h,u_h,iter_pf]=eulero_indietro(f,df,t_max,y_0,h)
% Risolve il problema di Cauchy
%
% y'=f(y,t)
% y(0)=y_0
%
% utilizzando il metodo di Eulero Indietro.
% Per ciascun istante temporale si calcola u_(n+1).
% Implementazione per il CASO SPECIFICO in cui la funzione del problema di
% Cauchy è LINEARE in y.
%
% Input:
% -> t_max: l'istante finale dell' intervallo temporale di soluzione
%                 (l'istante iniziale e' t_0=0)
% -> y_0: il dato iniziale del problema di Cauchy
% -> h: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% -> t_h: vettore degli istanti in cui si calcola la soluzione discreta
% -> u_h: la soluzione discreta calcolata nei nodi temporali t_h


t0 = 0;
t_h = t0:h:t_max;
% inizializzo il vettore che conterra' la soluzione discreta
N   = length(t_h);
u_h = zeros(1,N);
u_h(1) = y_0;
flineare = @(t) cos(t).*exp(-t/2);

for it=2:N
    u_old=u_h(it-1);
    t_new = t_h(it);
    u_h(it)=(u_old+h*flineare(t_new))/(1+0.5*h);
end