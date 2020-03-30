% Lab 7 Es 2

%1
lambda=-2; % coeff. problema modello
y0=1; % dato iniziale
tmax=10; % istante finale
t=linspace(0,10,100); 
h1=0.1; % passo di discretizzazione

f = @(t,y) lambda*y;
df= @(t,y) lambda;

f_esatta = @(t) y0*exp(lambda.*t); % soluzione esatta
y_esatta = f_esatta(t); % valori esatti

[tdisEA,ea]=eulero_avanti(f,t(end),y0,h1);
[tdisEI,ei]=eulero_indietro(f,df,t(end),y0,h1);
[tdisCN,cn]=crank_nicolson(f,df,t(end),y0,h1);

figure(1)
plot(tdisEA,ea,'bo-',tdisEI,ei,'ro-',tdisCN,cn,'go-',t,y_esatta,'k')
legend('Eulero avanti','Eulero indietro','Crank-Nicolson','F esatta')


%2
h2=0.9; % passo di discretizzazione

[tdisEA,ea]=eulero_avanti(f,t(end),y0,h2);
[tdisEI,ei]=eulero_indietro(f,df,t(end),y0,h2);
[tdisCN,cn]=crank_nicolson(f,df,t(end),y0,h2);

figure(2)
plot(tdisEA,ea,'bo-',tdisEI,ei,'ro-',tdisCN,cn,'go-',t,y_esatta,'k')
legend('Eulero avanti','Eulero indietro','Crank-Nicolson','F esatta')

%3
h3=1.1; % passo di discretizzazione

[tdisEA,ea]=eulero_avanti(f,t(end),y0,h3);
[tdisEI,ei]=eulero_indietro(f,df,t(end),y0,h3);
[tdisCN,cn]=crank_nicolson(f,df,t(end),y0,h3);

figure(3)
plot(tdisEA,ea,'bo-',tdisEI,ei,'ro-',tdisCN,cn,'go-',t,y_esatta,'k')
legend('Eulero avanti','Eulero indietro','Crank-Nicolson','F esatta')