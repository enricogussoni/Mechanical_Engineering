%Es 1 Lab 7

% Problema di Cauchy
f = @(t,y) cos(t).*exp(-t./2)-0.5*y;
df = @(t,y) -0.5*y;
y0=0;
t=linspace(0,10,1000);
f_es = @(t) sin(t).*exp(-t./2);

%1
% y_es=f_es(t);
% figure(1)
% plot(t,y_es)
% title('Soluzione esatta')

%3
h=0.5;
[ea,tea]=eulero_avanti(f,t(end),y0,h);
[ei,tei,inei]=eulero_indietro(f,df,t(end),y0,h);
[cn,tcn,incn]=crank_nicolson(f,df,t(end),y0,h);

figure(2)
plot(ea,tea,ei,tei,cn,tcn)
legend('Eulero avanti','Eulero indietro','Crank-Nicolson')