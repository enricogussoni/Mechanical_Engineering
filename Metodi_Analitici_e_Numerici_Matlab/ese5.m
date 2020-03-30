%% Esercizio 5

%% 1
fun=@(x) sin((1)./(1+x.^2));
a=-2*pi;
b=2*pi;
x_dis=linspace(a,b,1000);
f_dis=fun(x_dis);
figure
plot(x_dis,f_dis,'m')
n=10;
h=(b-a)/n;
x_nod=[a:h:b]; %11 nodi
f_nod=fun(x_nod);
grado=length(x_nod)-1;
P=polyfit(x_nod,f_nod,grado)
P_dis = polyval( P, x_dis );
hold on
plot(x_dis,P_dis,'g',x_nod,f_nod,'o')
xlabel('asse x')
ylabel('asse y')
title('Interpolazione f(x)')
legend('fun','interp','nodi')

pause

%% 2

err_dis = abs( P_dis - f_dis );
figure
plot(x_dis,err_dis)
xlabel('asse x')
ylabel('asse y')
title('Errore di interpolazione')
err_max=max(err_dis)

pause


%% 3

figure(3)
hold on
xlabel('asse x');
ylabel('asse y');
title('Interpolazione nodi Chebyshev');


figure(4)
hold on
xlabel('asse x')
ylabel('asse y')
title('Errore di interpolazione nodi Chebyshev')

N=[4 8 10];
for n=N
k=0:n;
t=-cos(pi*k/n);            % vettore contenete i nodi di interpolazione
x_nod=((b-a)/2)*t+(a+b)/2; % trasformazione affine
f_nod=fun(x_nod);
grado=n;                   % poniamo direttamente grado = n
P=polyfit(x_nod,f_nod,grado);
poly_dis=polyval(P,x_dis);
figure(3)
plot(x_dis,f_dis,'m',x_dis,poly_dis,'g',x_nod,f_nod,'o')
err_dis=abs(poly_dis-f_dis);
figure(4)
plot(x_dis,err_dis,'b')
fprintf('Errore con n=%d: %e\n',n,max(err_dis))
end
figure(3)
legend('fun','interp','nodi');

