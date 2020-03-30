% ESERCIZIO 3

clear all
close all
clc

% Grafico della funzione
fun = @(x) exp(-x.^2).*sin(x);
a = -2;
b =  3;
x_dis = linspace(a, b, 100);
f_dis = fun(x_dis);

figure(11)
plot(x_dis, f_dis, 'k')

%%

n = 30;
h = ( b - a ) / n;
x_nod = [ a : h : b ];

f_nod = fun(x_nod);

P_dis = interp1( x_nod, f_nod, x_dis );

hold on;
plot(x_nod,f_nod,'ro')

plot(x_dis,P_dis,'b*-')

legend('f(x)','nodi','interpolante composito lineare')
err   = max(abs(f_dis-P_dis))
%% PUNTO 3
H = [];
err_max = [];
for n = 2.^(2:7)
    h = ( b - a ) / n;
    H = [H h];
    x_nod = [ a : h : b ];
    f_nod = fun(x_nod);
    P_dis = interp1( x_nod, f_nod, x_dis );
    err_max = [ err_max; max( abs( P_dis - f_dis ) ) ];
end
figure(2)
semilogy(H,err_max,'ro-',H,H,'k--',H,H.^2,'k')
legend('errore','H','H^2')

