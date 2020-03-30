function [x, R, niter] = newtonsysGUS(F, JF, x0, tol, nmax)
% La funzione prende in ingresso la funzione F che restituisce la valutazione di F, la funzione JF
% che restituisce la valutazione dello Jacobiano della funzione F, il dato iniziale x0, il valore della
% tolleranza necessario per valutare il criterio d'arresto e il numero massimo di iterazioni.
% In uscita la funzione restituisce il vettore x contenente la soluzione del sistema non lineare, il
% valore del residuo valutato in corrispondenza della soluzione finale e il numero delle iterazioni
% effettuate.

niter=0;
x=x0;
R=tol+1; %per permettere il controllo nel 'while'

while (R>tol && niter<nmax)
    if (det(JF(x))==0)
        disp('Termine iterazione per annullamento jacobiana')
        break
    end
       delta=-JF(x) \ F(x);
       x=x+delta;
       R=norm(delta);
       niter=niter+1;
end

... disp comunicazioni su interruzione ciclo ecc...