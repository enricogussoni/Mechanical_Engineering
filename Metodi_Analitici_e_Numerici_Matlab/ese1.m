clc; clear all; close all;

fun = @(x) (x/2).*cos(x);
a = -2; b = 6;
x_dis = linspace( a, b,1000 );
f_dis = fun(x_dis);
figure(1)
plot( x_dis, f_dis, 'k' )
xlabel('x'); ylabel('y')
title('f(x) = x/2*cos(x)')
%%
n=2;                        % grado richiesto del polinomio
h=(b-a)/n;                  % passo tra i nodi
x_nod=[a:h:b];% vettore dei nodi di interpolazione
f_nod = fun(x_nod);         % valutazione della funzione nei nodi   

P = polyfit(x_nod,f_nod,n)
%
%x_dis=linspace(a,b,1000);

P_dis=polyval(P,x_dis);

hold on;
plot(x_nod,f_nod,'bo')
plot(x_dis,P_dis,'r')
title('f(x)=x/2*cos(x) e polinomio interpolatore di secondo grado')
legend('y = x/2 cos(x)','nodi equispaziati','pol 2')

%%
err_dis=abs( P_dis - f_dis );
figure(2)
plot(x_dis,err_dis);
grid
title('Errore di interpolazione')
err_max=max(err_dis)
%%

PP_dis = []; err_dis = []; err_max = [];
for n = [2 4 6]
    h = ( b - a ) / n;
    x_nod = [ a : h : b ];
    f_nod = fun(x_nod);
    P = polyfit( x_nod, f_nod, n );
    P_dis = polyval( P, x_dis );
    PP_dis = [ PP_dis; P_dis ];
    err_dis = [ err_dis; abs( P_dis - f_dis ) ];
    err_max = [ err_max; max( abs( P_dis - f_dis ) ) ];
end
figure
plot( x_dis, f_dis, 'k')
hold on
plot( x_dis, PP_dis )
legend('y=x/2 cos(x)','pol 2', 'pol 4', 'pol 6')
figure
plot(x_dis,err_dis)
grid
legend('err 2', 'err 4', 'err 6')