function [mat]=posvet(traccia,finder,vettore)

% traccia è il vettore dove vuoi avere i dati; finder è il vettore di cui
% 'vettore' è funzione e da cui vuoi ricavare le posizioni

for ii=1:length(traccia)
    fff=find(abs(finder-traccia(ii))<1e-5);
    mat(ii)=vettore(fff(end));
end
end
