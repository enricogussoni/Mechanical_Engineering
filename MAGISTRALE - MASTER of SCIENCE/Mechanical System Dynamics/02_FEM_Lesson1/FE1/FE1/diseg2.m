function diseg2(modo,fscala,incidenze,l,gamma,posiz,idb,xy);

%keyboard

% verifiche di consistenza dati in ingresso
[n_el,nc]=size(incidenze);

if length(posiz) ~= n_el
    sprintf('Errore: numero nodi nella matrice posiz diverso da numero elementi')
    return
end

n_gdl=length(modo);


hold on
% ciclo sugli elementi finiti
for k=1:n_el
% costruzione vettore spostamenti nodali dell'elemento nel s.d.r. globale
    for iri=1:6
        if incidenze(k,iri) <= n_gdl
        xkG(iri,1)=modo(incidenze(k,iri));
        else
        xkG(iri,1)=0.;
        end
    end
% applicazione fattore di scala
    xkG=fscala*xkG;
% proiezione da globale a locale vettore spostamenti nodali
    lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
              -sin(gamma(k)) cos(gamma(k)) 0.
	                0.      0.     1. ] ;
    Lambda = [ lambda     zeros(3)
              zeros(3)   lambda      ] ;
    xkL=Lambda*xkG;

% calcolo spostamenti assiali u e trasversali w mediante funzioni di forma
    csi=l(k)*[0:0.05:1];
    fu=zeros(6,length(csi));
    fu(1,:)=1-csi/l(k);
    fu(4,:)=csi/l(k);
    u=(fu'*xkL)';
    fw=zeros(6,length(csi));
    fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
    fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
    fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
    fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
    w=(fw'*xkL)';
% proiezione da s.d.r. locale a globale della deformata dell'elemento
    xyG=lambda(1:2,1:2)'*[u+csi;w];
    undef=lambda(1:2,1:2)'*[csi;zeros(1,length(csi))];
 % disegno traccia elemento deformato e indeformato
    plot(undef(1,:)+posiz(k,1),undef(2,:)+posiz(k,2),'b--')
    plot(xyG(1,:)+posiz(k,1),xyG(2,:)+posiz(k,2),'b')
end

% ciclo sui nodi
n_nodi=size(idb,1);
xkG=[];
for k=1:n_nodi
    for ixy=1:2
        if idb(k,ixy) <= n_gdl
        xkG(k,ixy)=modo(idb(k,ixy));
        else
        xkG(k,ixy)=0.;
        end
    end
end
xkG=fscala*xkG;
xyG=xkG+xy;
plot(xy(:,1),xy(:,2),'b.')
plot(xyG(:,1),xyG(:,2),'bo')

grid on
axis equal