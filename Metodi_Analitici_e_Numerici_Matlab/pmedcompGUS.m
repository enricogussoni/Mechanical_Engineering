function I = pmedcompGUS(a,b,N,fun)
% Implementa la formula di quadratura del punto medio
% composita su intervalli equispaziati. La funzione pmedcomp.m riceve in 
% ingresso gli estremi di integrazione a,b, il numero di sottointervalli N 
% in cui si vuole suddividere il dominio di integrazione e la funzione fun da
% integrare, e restituisce in uscita il valore approssimato dell'integrale.

H=(b-a)/N; % ampiezza sottointervalli
n=[a:H:b]; % nodi di interpolazione

% area retangolo = f(Mk)*H
% area integrale_approssimato = sum(area rettangolo)
% Mk = (x(k-1)+x(k))/2

m=(n(1:length(n)-1)+n(2:length(n)))./2; % punti medi dei sottointervalli
A=fun(m).*H; % aree dei singoli rettangoli
I=sum(A); % integrale approssimato