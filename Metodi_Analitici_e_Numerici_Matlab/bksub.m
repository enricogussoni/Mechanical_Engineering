function x=bksub(U,y)

% function [x] = bksub(U,y)
% Algoritmo di sostituzione all'indietro
%
% U: matrice quadrata triangolare superiore
% y: termine noto
% x: soluzione del sistema Ux=y

n=length(y);

if ((size(U,1) ~= n) || (size(U,2) ~= n))
  error('ERRORE in bksub: dimensioni incompatibili')
end

if (U~=triu(U))
  error('ERRORE in bksub: matrice non triangolare superiore')
end

if (prod(diag(U)) == 0)
  error('ERRORE in bksub: matrice singolare')
end

x = zeros(n,1);
x(n) = y(n) / U(n,n);
for i=n-1:-1:1
	x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
end
