function y=fwsub(L,b)

% function [y] = fwsub(L,b)
% Algoritmo di sostituzione in avanti
%
% L: matrice quadrata triangolare inferiore
% b: termine noto
% y: soluzione del sistema Ly=b

n=length(b);

if ((size(L,1) ~= n) || (size(L,2) ~= n))
  error('ERRORE in fwsub: dimensioni incompatibili')
end

if (L~=tril(L))
  error('ERRORE in fwsub: matrice non triangolare inferiore')
end

if (prod(diag(L)) == 0)
  error('ERRORE in fwsub: matrice singolare')
end

y = zeros(n,1);
y(1) = b(1) / L(1,1);
for i=2:n
	y(i) = (b(i) - L(i,1:i-1) * y(1:i-1)) / L(i,i);
end
