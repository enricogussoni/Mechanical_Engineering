function [x, niter, res, incr] = gradiente (A, b, x0, maxiter, toll)
%GRADIENTE Risolve un sistema lineare con il metodo del gradiente.
%   [X, NITER, RES, INCR] = GRADIENTE(A, B, X0, MAXITER, TOLL) 
%    risolve il sistema lineare A*X = B utilizzando il metodo del
%    gradiente, partendo da X0 e iterando fino a quando il residuo 
%    normalizzato non e' inferiore della tolleranza TOLL
%    o il numero massimo MAXITER di iterazioni e' superato.
%
%    Viene restituita la soluzione X, il numero di iterazioni NITER
%    effettuate, il vettore RES dei residui normalizzati e 
%    il vettore INCR degli incrementi relativi.

%% 1. Controllo che la matrice sia quadrata
[n, m] = size(A);
if (n ~= m) || (n ~= length(b))
  error('Dimensioni incompatibili')
end

%% 2. Algoritmo principale
% Inizializzo il numero di iterazioni, residuo ed errore
res  = zeros(maxiter, 1);
incr = zeros(maxiter, 1);
% Iterazioni
x = x0;
r = b - A*x;
for k = 1:maxiter
  % Salvo
  xold = x;
  % Passo k
  alpha = (r' * r) / (r' * A * r);
  x = xold + alpha * r;
  % Residuo
  r = b - A * x;
  res(k) = norm(r) / norm(b);
  % Errore
  incr(k) = norm(x - xold) / norm(xold);
  % Criterio di arresto
  if (res(k) < toll)
    break;
  end
end

% Iterazioni effettuate
niter = k;
res   = res(1:k);
incr  = incr(1:k);

%% 3. Info
if (niter < maxiter)
  fprintf('\n Convergenza raggiunta in %d iterazioni.\n', niter);
else
  fprintf('\n Convergenza NON raggiunta in %d iterazioni.\n', niter);
end

end % of function gradiente
