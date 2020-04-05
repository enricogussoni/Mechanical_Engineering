function [mat]=posvet(traccia,finder,vettore)
%traccia è il vettore dove vuoi avere i dati; finder è il vettore di cui
%'vettore' è funzione e da cui vuoi ricavare le posizioni%
for i=1:length(traccia)
    mat(i)=vettore(finder==traccia(i));
    %può non riconoscere uguali due numeri causa errore numerico, nel caso usare abs(finder-traccia(i))<1e-10
end
end
