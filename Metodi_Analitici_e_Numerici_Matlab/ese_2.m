
clc; clear; close all
%% punto 1
format short e
Iex = [1/30,  1/90, 1/182, 1/306];
c=1;
f =@(x) (1-x).*x.^(4*c);
a = 0; b = 1;
[x_nod,w]=gauleg(a,b,5);
I_n=sum(f(x_nod).*w)
I=1/30

figure; x=[a:0.01:b];
plot(x,f(x),'b-',x_nod,f(x_nod),'ro');

%% punto 2

f =@(x) (1-x).*x.^(4*c);

I=[];
for n=1:12
    [x_nod,w]=gauleg(a,b,n);
    I=[I sum(f(x_nod).*w)];
end
err = abs(Iex(c) - I);
err
figure; semilogy(err,'o');
xlabel('N')
ylabel('Errore')

%% punto 3 cambiare c e ripetere i punti 1 e 2
