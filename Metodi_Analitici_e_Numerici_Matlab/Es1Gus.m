% Lab 7 Es 1

% Problema di Caucy
f_t_y = @(t,y) cos(t).*exp((-t)./2)-0.5.*y;
df_y = @(t,y) -1/2;
y0=0; % Condizione al contorno
t=linspace(0,10,1000); % Intervallo

%1
sol_esatta = @(t) sin(t).*exp(-t./2);
y_esatta = sol_esatta(t);
figure(1)
plot(t,y_esatta)
xlabel('Tempo')
ylabel('Soluzione esatta')
title('Soluzione del PdC')
hold on

%3
h=0.5; % passo di discretizzazione

[t_hea,EA]=eulero_avanti(f_t_y,t(end),y0,h);
[t_hei,EI]=eulero_indietro(f_t_y,df_y,t(end),y0,h);
[t_hcn,CN]=crank_nicolson(f_t_y,df_y,t(end),y0,h);

figure(2)
plot(t_hea,EA,'bo-',t_hei,EI,'ro-',t_hcn,CN,'go-',t,y_esatta,'k')
legend('Eulero avanti','Eulero indietro','Crank-Nicolson')

%4
h=[0.4, 0.2, 0.1, 0.05, 0.025, 0.0125];
n=length(h);
% ea=zeros(1,n);
% ei=zeros(1,n);
% cn=zeros(1,n);

err_EA=zeros(1,n);
err_EI=zeros(1,n);
err_CN=zeros(1,n);

for i=1:n
    [tea,ea]=eulero_avanti(f_t_y,t(end),y0,h(i));
    err_EA = max(abs(ea - sol_esatta(tea)));
    
    [tei,ei]=eulero_indietro(f_t_y,df_y,t(end),y0,h(i));
    err_EI = max(abs(ei - sol_esatta(tei)));
    
    [tcn,cn]=crank_nicolson(f_t_y,df_y,t(end),y0,h(i));
    err_CN = max(abs(cn - sol_esatta(tcn)));
end

figure(3)
loglog(h,err_EA,h,err_EI,h,err_CN,h,h,h,h.^2)
xlabel('Passo di discretizzazione')
ylabel('Errore')
legend('Eulero Avanti','Eulero Indietro','Crank Nicolson','h','h^2')
