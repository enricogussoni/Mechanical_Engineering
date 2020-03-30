clear all; close all
N = 10;
n = 1:N;

%************ DST *************

f = @(x) sin(x) + 0.5*sin(2*x) + 0.2*sin(5*x);


% calcolo coefficienti bk tramite dst
xs = pi*n/(N+1);
bk = 2/(N+1)*dst(f(xs));

% plotto i coefficienti bk
figure(1)
plot(n,bk,'o-')

%ricostruisco la funzione
xdis = linspace(-pi,pi,1000);
fric=0;
for i=n
    fric=fric+bk(i)*sin(i*xdis);
end

%confronto
figure(2)
plot(xdis,f(xdis),'r','Linewidth',4); hold on; grid
plot(xdis,fric,'b','Linewidth',2)
legend('f(x)','fric')

%%
close all
%************ DCT *************

g = @(x) cos(x)+3*cos(2*x)+2;

% calcolo coefficienti ak tramite dct
xc=  pi*(2*n - 1)/(2*N);
vec_dct=dct(g(xc));

a0=2*sqrt(1/N)*vec_dct(1);
ak=sqrt(2/N)*vec_dct(2:end);

ak=[a0 ak];

% plotto i coefficienti ak
figure(1)
plot(n-1,ak,'o-')

%ricostruisco la funzione
xdis = linspace(-pi,pi,1000);
gric = ak(1)/2;
for i=2:N
    gric=gric+ak(i)*cos((i-1)*xdis);
end

%confronto
figure(2)
plot(xdis,g(xdis),'r','Linewidth',4); hold on; grid
plot(xdis,gric,'b','Linewidth',2)
legend('g(x)','gric')

%%
close all
%Funzione da scomporre in pari + dispari
h = @(x) sin(x) + 0.5*cos(2*x);

hp= @(x) 0.5*(h(x) + h(-x));
hd= @(x) 0.5*(h(x) - h(-x));

% calcolo coefficienti bk tramite dst
xs = pi*n/(N+1);
bk = 2/(N+1)*dst(hd(xs));

% calcolo coefficienti ak tramite dct
xc=  pi*(2*n - 1)/(2*N);
vec_dct=dct(hp(xc));

a0=2*sqrt(1/N)*vec_dct(1);
ak=sqrt(2/N)*vec_dct(2:end);

ak=[a0 ak];

%plot
figure(1)
plot(n,bk,'o-')
legend('bk')
figure(2)
plot(n-1,ak,'o-')
legend('ak')
%ricostruisco la funzione
xdis = linspace(-pi,pi,1000);
hric = ak(1)/2;
for i=2:N
    hric=hric+ak(i)*cos((i-1)*xdis);
end
for i=n
    hric=hric+bk(i)*sin(i*xdis);
end

figure(3)
plot(xdis,h(xdis),'r','Linewidth',4); hold on; grid
plot(xdis,hric,'b','Linewidth',2)
legend('h(x)','hric','Location','se')

%%
%Funzione qualsiasi da prolungare pari o dispari:
close all
clear all; close all
N = 10;
n = 1:N;

%Prolumgamento dispari = utilizzo dst
wd = @(x) x;
% calcolo coefficienti bk tramite dst
xs = pi*n/(N+1);
bk = 2/(N+1)*dst(wd(xs));

% plotto i coefficienti bk
figure(1)
plot(n,bk,'o-')
legend('bk')

%ricostruisco la funzione
xdis = linspace(-pi,pi,1000);
wricd=0;
figure(2)
for i=n
    wricd=wricd+bk(i)*sin(i*xdis);
end
%confronto
plot(xdis,wd(xdis),'r','Linewidth',4); hold on; grid
plot(xdis,wricd,'b','Linewidth',2)
legend('w(x)','wricd','Location','se')

%% Prolungo pari = utilizzo dct
close all
wp = @(x) abs(x);


% calcolo coefficienti ak tramite dct
xc = pi*(2*n - 1)/(2*N);
vec_dct = dct(wp(xc));

a0 = 2*sqrt(1/N)*vec_dct(1);
ak = sqrt(2/N)*vec_dct(2:end);

ak=[a0 ak];

% plotto i coefficienti ak
figure(3)
plot(n-1,ak,'o-')
legend('ak')

%ricostruisco la funzione
xdis = linspace(-pi,pi,1000);
wricp = ak(1)/2;
figure(4)
for i=2:N
    wricp=wricp+ak(i)*cos((i-1)*xdis);
end

%confronto
plot(xdis,wp(xdis),'r','Linewidth',4); hold on; grid
plot(xdis,wricp,'b','Linewidth',2)
legend('w(x)','wricp','Location','se')
