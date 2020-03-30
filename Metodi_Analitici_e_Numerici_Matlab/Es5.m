% Lab 5 Es 5

f = @(x) sin(1./(1+x.^2)); % funzione
x = linspace(-2*pi, 2*pi, 100); % punti nell'intervallo
y = f(x); % valori di f nei punti x

%1 polinomio di Lagrange
n=10; % grado del polinomio
a=x(1); % estremo sinistro dell'intervallo
b=x(length(x)); % estremo destro dell'intervallo
h=(b-a)/n; % ampiezza sottointervalli
nodi=[a:h:b]; % nodi di interpolazione (11)
p=polyfit(nodi,f(nodi),length(nodi)-1); % coefficienti del polinomio di Lagrange
P=polyval(p,nodi); % polinomio inerpolante

% stampa risultato
disp('I coefficienti del polinomio interpolante sono')
p

% grafico di controllo
figure(1)
plot(x,f(x),nodi,P,'r')
xlabel('x')
ylabel('funzione e interpolante')

%2 valutazione errore
err=max(abs((f(nodi)-P))); % errore in norma infinito
                           % valutato nei nodi
                           
%3 %SBAGLIATO (vedi correzioni punto 1)
% figure(2)
% plot(x,y)
% hold on
c=[4 8 10];
k=1;
 
while (k<4)
    n=c(k);
    I=0:n; % nodi base equispaziati
    nodic=-cos(pi*I./n); % nodi C-G-L di base
    nodic_def=(a+b)/2 + ((b-a)/2).*nodic; % nodi C-G-L traslati
    
    nnodi=11; % numero nodi di interpolazione
    p=polyfit(x,y,nnodi-1); % coefficienti del polinomio di Lagrange
    P=polyval(p,nodic_def); % polinomio inerpolante
    
    err(k)=max(abs((f(nodic_def)-P)));
    
    figure(k)
    plot(nodic_def,P)
    
    k=k+1;
end

% hold off

figure(3)
plot(c,err)
title('Riduzione errore')
xlabel('Grado del polinomio interpolante')
ylabel('Errore')