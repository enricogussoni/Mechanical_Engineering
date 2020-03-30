function [T, niter] = iterqr (A, maxiter, toll)
%ITERQR Iterazioni QR.
%   [T, ERR] = ITERQR(A, MAXITER, TOLL)
%    Calcola la matrice T con le iterazioni QR. Da quest'ultima
%    si possono poi leggere gli autovalori approssimati di A sulla
%    diagonale.
%    L'algoritmo si ferma quando la norma 2 della sottodiagonale
%    e' inferiore a TOLL, oppure se il numero massimo MAXITER
%    di iterazioni e' stato raggiunto.

%% 1. Algoritmo principale
% Inizializzo
T = A;
% Iterazioni
for k = 1:maxiter
    % Fattorizzazione QR
    [Q, R] = qr(T);
    T = R * Q;
    % Norma della sottodiagonale
    if (norm(diag(T, -1)) < toll)
        break;
    end
end
niter = k;

if (niter < maxiter)
    fprintf('\n Convergenza raggiunta in %d iterazioni.\n', niter);
else
    fprintf('\n Convergenza NON raggiunta in %d iterazioni.\n', niter);
end

end % of function iterqr
