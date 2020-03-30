% Lab 8 Es 1 (2)

N=10;
n=1:N;

%f
x = linspace(-pi,pi,100);
f = @(x) sin(x) + 0.5*sin(2*x) + 0.2*sin(5*x);

xf=f(pi*n./(N+1));
yf=dst(xf);
bf=yf.*2/(N+1);

figure(1)
plot(n,bf)
title('Coeff. di Fourier di f')

f_esatta = f(x);
f_ricostruita = 0;

for c=n
    f_ricostruita = f_ricostruita + bf(c)*sin(c*x); 
end

figure(12)
plot(x,f_esatta,x,f_ricostruita)
title('f')
legend('f esatta','f ricostruita')

%g
g = @(x) cos(x) + 3*cos(2*x) +2;

xg=g(pi*(2*n-1)/(2*N));
yg=dct(xg);
ag=zeros(1,length(yg));
ag(1)=yg(1)*2/sqrt(N);
ag(2:end)=yg(2:end).*sqrt(2/N);

figure(2)
plot(n-1,ag)
title('Coeff. di Fourier di g')

g_esatta = g(x);
g_ricostruita = ag(1)/2;

for c=2:N
    g_ricostruita = g_ricostruita + ag(c)*cos((c-1)*x); 
end

figure(22)
plot(x,g_esatta,x,g_ricostruita)
title('g')
legend('g esatta','g ricostruita')

%h
h = @(x) sin(x) + 0.5*cos(2*x);
ph = @(x)(h(x)+h((-1)*x))/2;
dh = @(x) (h(x)-h((-1)*x))/2;

xdh=dh(pi*n./(N+1));
ydh=dst(xdh);
bdh=ydh.*2/(N+1);

xph=ph(pi*(2*n-1)/(2*N));
yph=dct(xph);
aph=zeros(1,length(yph));
aph(1)=yph(1)*2/sqrt(N);
aph(2:end)=yph(2:end).*sqrt(2/N);

figure(3)
plot(n,bdh,'o-',n-1,aph)
legend('Coeff. parte dispari','Coeff. parte pari');
title('Coeff. di Fourier di h');

h_esatta = h(x);

dh_esatta = dh(x);
dh_ricostruita = 0;

for c=n
    dh_ricostruita = dh_ricostruita + bdh(c)*sin(c*x); 
end

ph_esatta = ph(x);
ph_ricostruita = aph(1)/2;

for c=2:N
    ph_ricostruita = ph_ricostruita + aph(c)*cos((c-1)*x); 
end

h_ricostruita = dh_ricostruita + ph_ricostruita;

figure(32)
plot(x,h_esatta,x,h_ricostruita);
title('h')
legend('h esatta','h ricostruita')

%w
%w = @(x) x;
wd = @(x) x;
wp = @(x) abs(x);

xwd=wd(pi*n./(N+1));
ywd=dst(xwd);
bwd=ywd.*2/(N+1);

xwp=wp(pi*(2*n-1)/(2*N));
ywp=dct(xwp);
awp=zeros(1,length(ywp));
awp(1)=ywp(1)*2/sqrt(N);
awp(2:end)=ywp(2:end).*sqrt(2/N);

figure(4)
plot(n,bwd,'o-',n-1,awp)
legend('Coeff. parte dispari','Coeff. parte pari');
title('Coeff. di Fourier di w');

w_esatta = x;

wd_esatta = wd(x);
wd_ricostruita = 0;

for c=n
    wd_ricostruita = wd_ricostruita + bwd(c)*sin(c*x); 
end

wp_esatta = wp(x);
wp_ricostruita = awp(1)/2;

for c=2:N
    wp_ricostruita = wp_ricostruita + awp(c)*cos((c-1)*x); 
end

w_ricostruita = wd_ricostruita + wp_ricostruita;

figure(42)
plot(x,w_esatta,x,w_ricostruita);
title('w')
legend('w esatta','w ricostruita') %???

