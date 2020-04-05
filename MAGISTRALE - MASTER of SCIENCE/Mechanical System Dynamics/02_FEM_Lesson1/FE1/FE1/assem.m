function [M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb)

% verifiche di consistenza dati in ingresso
[n_el,nc]=size(incidenze);
if nc ~= 6
    disp('Errore: numero colonne matrici incidenze diverso da 6')
    return
end

if length(l) ~= n_el
    sprintf('Errore: numero parametri l diverso da numero elementi')
    return
end
if length(m) ~= n_el
    sprintf('Errore: numero parametri m diverso da numero elementi')
    return
end
if length(EA) ~= n_el
    sprintf('Errore: numero parametri EA diverso da numero elementi')
    return
end
if length(EJ) ~= n_el
    sprintf('Errore: numero parametri EJ diverso da numero elementi')
    return
end
if length(gamma) ~= n_el
    sprintf('Errore: numero parametri gamma diverso da numero elementi')
    return
end

if min(min(incidenze)) ~= 1
    sprintf('Errore: la numerazione dei gradi di libertà non incomincia da 1')
    return
end

% limite sui gdl totali struttura
n_gdl=max(max(idb));
if n_gdl > 100
    sprintf('Errore: numero g.d.l. > 100 non permesso')
    return
end

% sprintf('Numero elementi , %d',n_el)
% sprintf('Numero g.d.l.   , %d',n_gdl)
% pause(1)

% assemblaggio marici M e K
M=zeros(n_gdl,n_gdl);
K=zeros(n_gdl,n_gdl);
for k=1:n_el
    [mG,kG] = el_tra(l(k),m(k),EA(k),EJ(k),gamma(k));
    for iri=1:6
        for ico=1:6
            i1=incidenze(k,iri);
            i2=incidenze(k,ico);
            M(i1,i2)=M(i1,i2)+mG(iri,ico);
            K(i1,i2)=K(i1,i2)+kG(iri,ico);
        end
    end
end