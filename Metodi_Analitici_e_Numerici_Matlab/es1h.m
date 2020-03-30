%%
clear
clc

h = @(x) sin(x) + 0.5*cos(2*x);
N=10;
n=1:N;
x0=-pi;
xN=pi;

xnodid = n*pi/(N+1);

xnodip = pi*(2*n-1)/(2*N);

% hd = @(x) sin(x);
% hp = @(x) 0.5*cos(2*x);
hd = @(x) 0.5*(h(x)-h(-x));
hp = @(x) 0.5*(h(x)+h(-x));

xd=hd(xnodid);
xp=hp(xnodip);

yd=dst(xd);
yp=dct(xp);

a0=(2/sqrt(N))*yp(1);
a=[a0 sqrt(2/N)*yp(2:end)];
b=(2/(N+1))*yd;

figure(2)
title('coefficienti')
plot(1:N,b,1:N,a)
legend('dispari','pari')

%%

x=linspace(x0,xN,1000);
hric=a0/2;
for c=2:N;
    hric=hric+a(c)*cos((c-1)*x);
end
for c=n
    hric=hric+b(c)*sin(c*x);
end

figure(3)
plot(x,h(x),x,hric,'o-')
