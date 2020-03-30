function [p, c] = stimap (xvect, toll)
%
% [p,c]=stimap(xvect,toll)
%
% Stima ordine e fattore di abbattimento di un metodo 
% iterativo per il calcolo degli zeri di una funzione 
% utilizzando le seguenti formule :
%          
%        | x_(k+1) - x_k |
%     ln ------------------
%        | x_k - x_(k-1) |            | x_(k+1) - x_k |
% p = --------------------------  c = ---------------------
%        | x_k - x_(k-1) |           | x_(k) - x_(k-1) |^p  
%     ln ----------------------
%        | x_(k-1) - x_(k-2) |
%
% Parametri di ingresso:
%
% xvect         Vettore contenente tutte le iterate
% toll          Soglia sotto la quale due iterate successive 
%               sono considerate troppo vicine
%
%          
% Parametri di uscita:
%
% p          Vettore contenente tutte le stime dell'ordine calcolate
% c          Vettore contenente tutte le stime del 
%            fattore di abbattimento dell'errore
  if (nargin < 2)
    toll = eps;
  end
  
  %% 1. xvect deve avere almeno 4 valori.
  if (length(xvect) < 4)
    error('xvect ha meno di 4 valori, stima non possibile!');
  end
  
  %% 2. Stima dell'errore
  err = abs( xvect(2:end) - xvect(1:end-1) );
  % Se ci sono valori molto piccoli, la stima e' mal condizionata
  if ( find(err < toll) )
    error('xvect contiene valori molto vicini nella tolleranza fissata!');
  end
  % Se per un metodo numerico iterativo vale asintoticamente
  % e_k+1 = c * (e_k)^p , si dice che il metodo ha ordine di convergenza p e fattore di convergenza c.
  % prendendo il logaritmo:
  %   log(e_{k+1}) = log(c) + p * log(e_{k})
  % Facciamo una stima ai minimi quadrati con polyfit.
  
  % log(e_{k+1})
  err_kp1 = log(err(2:end));
  % log(e_{k})
  err_k   = log(err(1:end-1));
  % Stima
  coeff = polyfit(err_k, err_kp1, 1);
  p = coeff(1);
  c = exp(coeff(2));
  % Stampo il risultato
  fprintf(' Ordine stimato       : %12.8f \n', p);
  fprintf(' Fattore di riduzione : %12.8f \n', c);
  
end




