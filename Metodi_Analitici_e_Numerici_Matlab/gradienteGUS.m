function [x, niter, res, incr] = gradienteGUS (A, b, x0, maxiter, toll)

% Lab 3 Es 1.2
% La funzione deve prendere in ingresso, oltre alla
% matrice A e il vettore b, il dato iniziale x0, il numero massimo di 
%iterazioni kmax e la tolleranza per il controllo sul residuo normalizzato.
% In uscita la funzione deve fornire la soluzione x0 del sistema lineare, 
% il numero di iterazioni effettuate niter e i vettori res e incr contenenti 
% rispettivamente i residui normalizzati e gli incrementi relativi calcolati 
% ad ogni iterazione.

y=zeros(maxiter,1); %allocare memoria serve anche a cautelarsi da loop
r=zeros(maxiter,1);
a=zeros(maxiter,1);
res=zeros(maxiter,1); %CORRETTO
incr=zeros(maxiter,1); %"
n=length(x0);
y=ones(n,maxiter);
y(1:16,1)=x0;

for k=1:maxiter %stop di sicurezza al maxiter° ciclo 
    r(k)=b-A*y(k);
    a(k)=(r(k)'*r(k))/(r(k)'*A*r(k));
    y(k+1)=y(k)+a(k)*r(k);
    res(k)=norm(r(k))/norm(b);
    incr(k)=norm(y(k+1)-x(k))/norm(x(k));
    %il ciclo, oltre che non funzionare, è inutile
    %perché non serve memorizzare ogni r e ogni x
    %ma solo quelli finali
    
    % controllo sull'errore
    if res(k)<toll
       niter=k;     % per finire il ciclo si può usare "break"
       k=maxiter;   % k non "muore" alla fine del ciclo ma rimane in memoria
    end
end

x=y(niter,1); %per usare solo i valori diversi da 0 di y
res=res(1:niter);
incr=incr(1:niter);