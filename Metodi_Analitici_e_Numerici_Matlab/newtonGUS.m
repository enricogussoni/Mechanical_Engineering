function [xvect,it]=newtonGUS(x0,nmax,toll,fun,dfun,mol)

% La funzione prende in ingresso il dato iniziale x(0), il numero massimo di iterazioni, il valore
% della tolleranza epxylon necessario per il criterio d'arresto, la funzione f di cui si stanno ricercando gli
% zeri, la sua derivata f' e per ultimo il valore della molteplicità della radice.
% In uscita la funzione restituisce il vettore xvect contenente tutte le iterate calcolate (l'ultima
% componente sarà perciò la soluzione) e il numero di iterazioni effettuate.

xvect=zeros(nmax,1);
xvect(1)=x0;

for c=1:nmax  %Inserire controllo su df diverso da 0
    xvect(c+1)=xvect(c)-mol*(fun(xvect(c))/dfun(xvect(c))); 
    it=c;
    if (xvect(c+1)-xvect(c)<toll) %Criterio di arresto
        break
    end
end

% Un ciclo while permette di velocizzare l'operazione (vedi soluzione)

xvect=xvect(1:it);