clc; clear all; close all

fun = @(x) 1./(1+x.^2);

n=15;
a=-5;
b=5;
h=(b-a)/n;
x_nod=[a:h:b];


f_nod = fun(x_nod);

grado=(length(x_nod))-1; % verifica del grado del pol. interpolatore
P=polyfit(x_nod,f_nod,grado)


x_dis=linspace(a,b,1000);
poly_dis = polyval( P, x_dis );

f_dis=fun(x_dis);
figure
plot(x_dis,f_dis,'m',x_dis,poly_dis,'g')
xlabel('asse x')
ylabel('asse y')
title('Interpolazione nodi equispaziati')
legend('fun','interp')

err_dis = abs( poly_dis - f_dis );
figure
xlabel('asse x')
ylabel('asse y')
title('Errore di interpolazione nodi equispaziati')
plot(x_dis,err_dis,'m')

%%
k = [0:n];
t = -cos(pi*k/n);             % vettore contenete i nodi di interpolazione
x_nod = (a+b)/2 + ((b-a)/2)*t;

f_nod=fun(x_nod);
grado=length(x_nod)-1;
P=polyfit(x_nod,f_nod,grado);
poly_dis = polyval( P, x_dis );
figure
plot(x_dis,f_dis,'m',x_dis,poly_dis,'g')
xlabel('asse x')
ylabel('asse y')
title('Interpolazione nodi di Chebyshev')
legend('fun','interp')

err_dis = abs( poly_dis - f_dis );
figure
plot(x_dis,err_dis,'b')
xlabel('asse x')
ylabel('asse y')
title('Errore di interpolazione nodi di Chebyshev')
    
err_dis = abs( poly_dis - f_dis );
err_max=max(err_dis)
