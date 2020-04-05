function [file_i,xy,nnod,sizew,idb,ngdl,incidenze,l,gamma,m,EA,EJ,posiz,nbeam]=loadstructure

% Caricamento file input (Si può modificare con uigetfile)
disp(' ');
file_i=input(' Nome file *.inp (senza est.) di input della struttura = ','s') ;
disp(' ')
% Verifica esistenza file
if exist([file_i '.inp'])~=2
    disp(['Il file ',file_i,'.inp non esiste' ])
    file_i=[];
    return
end

%nprova=input(' Nome della prova = ','s');
nprova=file_i;

% Apro il file.
eval(['fid_i=fopen(''',file_i,'.inp'',''r'');']);

% Lettura card nodi
findcard(fid_i,'*NODES')
% line=scom(fid_i)
% finewhile=findstr(line,'*ENDNODES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDNODES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i %f %f')';
        iconta=iconta+1;
        if iconta ~=tmp(1)
        disp('Errore: nodi non numerati in ordine progressivo')
        break  
        end
        ivinc(iconta,:)=tmp(2:4);
        xy(iconta,:)=tmp(5:6);
    end
end
nnod=iconta;
disp(['Numero nodi struttura ',int2str(nnod)])
sizew=sqrt((max(xy(:,1))-min(xy(:,1)))^2+(max(xy(:,2))-min(xy(:,2)))^2);
% fine lettura card nodi

% Costruzione matrice IDB di numerazione gdl nodi
% 1) g.d.l. liberi
igdl=0;
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 0
            igdl=igdl+1;
            idb(i,j)=igdl;
        end
    end
end
ngdl=igdl;
disp(['Numero g.d.l. struttura ',int2str(ngdl)])
% 2) g.d.l. vincolati
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 1
            igdl=igdl+1;
            idb(i,j)=igdl;
        end
    end
end


% Lettura card travi
findcard(fid_i,'*BEAMS')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDBEAMS'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %f %f %f')';
        iconta=iconta+1;
        incid(iconta,:)=tmp(2:3);
        m(1,iconta)=tmp(4);
        EA(1,iconta)=tmp(5);
        EJ(1,iconta)=tmp(6);
        l(1,iconta)=sqrt((xy(incid(iconta,2),1)-xy(incid(iconta,1),1))^2+(xy(incid(iconta,2),2)-xy(incid(iconta,1),2))^2);
        gamma(1,iconta)=atan2(xy(incid(iconta,2),2)-xy(incid(iconta,1),2),xy(incid(iconta,2),1)-xy(incid(iconta,1),1));
        incidenze(iconta,:)=[idb(incid(iconta,1),:) idb(incid(iconta,2),:)];
        posiz(iconta,:)=xy(incid(iconta,1),:);
    end
end
nbeam=iconta;
disp(['Numero elementi trave ',int2str(nbeam)])
% fine lettura card nodi

% Lettura card "Damping"
% findcard(fid_i,'*DAMPING')
% line=scom(fid_i);
% tmp=sscanf(line,'%f %f')';
% alfa=tmp(1);
% beta=tmp(2);
% fine lettura card "Damping"

% Chiudo il file della struttura.
fclose(fid_i) ;