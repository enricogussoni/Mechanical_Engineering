%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Corso di Metodi di Calcolo delle Strutture       %
%               Analisi statica di telai 2D              %
%                      13/01/2016                        %
%                                     Giuseppe COCCHETTI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
clc
close all


db=1000; %[mm]

%Coordinate dei nodi
dXY=[   0, db;
        0,  0;
       db,  0;
     2*db, db];

nNodi=size(dXY,1);  %Numero totale dei nodi
  


%Aste e connessioni (tabelle delle incidenze)
nIncNodi=[1,2;
          2,3;
          3,4];

%nIncGdl=[ 1, 2, 3, 4, 5, 6;
%          4, 5, 6, 7, 8, 9;
%         13, 8, 9,10,11,12];
nIncGdl=[nIncNodi(:,1)*3-2, nIncNodi(:,1)*3-1, nIncNodi(:,1)*3, nIncNodi(:,2)*3-2, nIncNodi(:,2)*3-1, nIncNodi(:,2)*3];
nIncGdl(3,1)=13;

nAste=size(nIncNodi,1);

nGdlTot=max(max(nIncGdl));   %Numero totale dei gradi di liberta' (gdl)




%Parametri meccanici (rigidezza assiale EA e rigidezza flessionale EI)
%    [   E  * A ,        E * I]
dEAI=[206000*(10*30),   206000*(10*30^3/12);
      206000*(10*30),   206000*(10*30^3/12);
      206000*(10*30),   206000*(10*30^3/12)];



%Gdl da vincolare 
nUv=[1,4,10,11]; 

%Gdl da non vincolare
%nUl=[2,3,5,6,7,8,9,12,13]; 
nUl=[2,3,5,6,7,8,9,12,13];
% nUl(nUv)=[];


%Vettore termini noti
dF=zeros([nGdlTot,1]);
dF(13,1)=dF(13,1)-5*10;

%Carichi distribuiti (sistema di riferimento locale) [p, q, m]
dq=-500/db;
dpqm=[0, 2*dq,0;
      0,    0,0;
      0,    0,0];

%Cedimenti vincolari assegnati
dUv=[0,0,0,0]';

%%%%%%%%%%%%%%%%%%%%%%% Termine assegnazione dati %%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Elaborazione %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Costruzione della matrice di rigidezza globale (K)
dK=zeros([nGdlTot,nGdlTot]);
for ne=1:nAste
  nA=nIncNodi(ne,1);  %nA=numero del primo nodo dell'asta
  nB=nIncNodi(ne,2);  %nB=numero del secondo nodo dell'asta

  dXYA=dXY(nA,:);  %dXYA=coordinate del primo nodo dell'asta
  dXYB=dXY(nB,:);  %dXYB=coordinate del secondo nodo dell'asta

  dL=norm(dXYB-dXYA);  %dL=lunghezza dell'asta

  dEAn=dEAI(ne,1); %Rigidezza assiale dell'asta analizzata (numero "ne")
  dEIn=dEAI(ne,2); %Rigidezza flessionale dell'asta analizzata (numero "ne")

  %Matrice di rigidezza della trave analizzata
  dk=[ dEAn/dL,             0,            0, -dEAn/dL,             0,            0;
             0,  12*dEIn/dL^3,  6*dEIn/dL^2,        0, -12*dEIn/dL^3,  6*dEIn/dL^2;
             0,   6*dEIn/dL^2,    4*dEIn/dL,        0,  -6*dEIn/dL^2,    2*dEIn/dL;
      -dEAn/dL,             0,            0,  dEAn/dL,             0,            0;
             0, -12*dEIn/dL^3, -6*dEIn/dL^2,        0,  12*dEIn/dL^3, -6*dEIn/dL^2;
             0,   6*dEIn/dL^2,    2*dEIn/dL,        0,  -6*dEIn/dL^2,    4*dEIn/dL]; 
  
  df=zeros([6,1]);
  dp=dpqm(ne,1); 
  dq=dpqm(ne,2); 
  dm=dpqm(ne,3); 
  df=[dp*dL/2, dq*dL/2-dm, dq*dL^2/12, dp*dL/2, dq*dL/2+dm, -dq*dL^2/12]';
  
  %Matrice di rotazione per l'elemento finito di trave
  de=(dXYB-dXYA)/dL;  %Versore assiale
  dCosAlpha=de(1);
  dSinAlpha=de(2);
  dRt=[ dCosAlpha, dSinAlpha, 0;
       -dSinAlpha, dCosAlpha, 0;
                0,         0, 1];

  dQ=[         dRt, zeros([3,3]);
      zeros([3,3]),         dRt];
              
  dk=dQ'*dk*dQ;
  df=dQ'*df;
  
  %Assemblaggio della matrice dk e del vettore df
  nv=nIncGdl(ne,:); %Vettore dei gdl relativi ai nodi dell'asta analizzata
  dK(nv,nv)=dK(nv,nv)+dk;
  dF(nv,1)=dF(nv,1)+df;
end



%%%%%%%%%%%%%%% Risoluzione del sistema lineare %%%%%%%%%%%%%%%%%%%
%Partizione matrici e vettori con riferimento ai gdl liberi e a quelli da vinvolare
dKll=dK(nUl,nUl);
dKlv=dK(nUl,nUv);
dKvl=dK(nUv,nUl);
dKvv=dK(nUv,nUv);

dFl=dF(nUl,1);
dFv=dF(nUv,1);

dUl=dKll\(dFl-dKlv*dUv);
dSv=dKvl*dUl+dKvv*dUv-dFv;

% Riposizionamento di U e S
dU=zeros([nGdlTot,1]);
dU(nUl,1)=dUl;
dU(nUv,1)=dUv;
dU

dS=zeros([nGdlTot,1]);
dS(nUv,1)=dSv;
dS



%%%%%%%%%%%%%%% Determinazione delle Azioni interne %%%%%%%%%%%%%%%%%%%
dNTM=zeros([6,nAste]);
for ne=1:nAste
  nA=nIncNodi(ne,1);  %nA=numero del primo nodo dell'asta
  nB=nIncNodi(ne,2);  %nB=numero del secondo nodo dell'asta

  dXYA=dXY(nA,:);  %dXYA=coordinate del primo nodo dell'asta
  dXYB=dXY(nB,:);  %dXYB=coordinate del secondo nodo dell'asta

  dL=norm(dXYB-dXYA);  %dL=lunghezza dell'asta

  
  dEAn=dEAI(ne,1); %Rigidezza assiale dell'asta analizzata (numero "ne")
  dEIn=dEAI(ne,2); %Rigidezza flessionale dell'asta analizzata (numero "ne")

  %Matrice di rigidezza della trave analizzata
  dk=[ dEAn/dL,             0,            0, -dEAn/dL,             0,            0;
             0,  12*dEIn/dL^3,  6*dEIn/dL^2,        0, -12*dEIn/dL^3,  6*dEIn/dL^2;
             0,   6*dEIn/dL^2,    4*dEIn/dL,        0,  -6*dEIn/dL^2,    2*dEIn/dL;
      -dEAn/dL,             0,            0,  dEAn/dL,             0,            0;
             0, -12*dEIn/dL^3, -6*dEIn/dL^2,        0,  12*dEIn/dL^3, -6*dEIn/dL^2;
             0,   6*dEIn/dL^2,    2*dEIn/dL,        0,  -6*dEIn/dL^2,    4*dEIn/dL]; 
  
  df=zeros([6,1]);
  dp=dpqm(ne,1); 
  dq=dpqm(ne,2); 
  dm=dpqm(ne,3); 
  df=[dp*dL/2, dq*dL/2-dm, dq*dL^2/12, dp*dL/2, dq*dL/2+dm, -dq*dL^2/12]';

  
  %Matrice di rotazione per l'elemento finito di trave
  de=(dXYB-dXYA)/dL;  %Versore assiale
  dCosAlpha=de(1);
  dSinAlpha=de(2);
  dRt=[ dCosAlpha, dSinAlpha, 0;
       -dSinAlpha, dCosAlpha, 0;
                0,         0, 1];

  dQ=[         dRt, zeros([3,3]);
      zeros([3,3]),         dRt];
  
  
  nv=nIncGdl(ne,:); %Vettore dei gdl relativi ai nodi dell'asta analizzata
  dUab=dU(nv,1);  %Spostamenti relativi ai nodi dell'asta analizzata (sistema globale)
  dUab=dQ*dUab; %Spostamenti relativi ai nodi dell'asta analizzata (sistema locale)
  dNTM(:,ne)=dk*dUab-df; %Azioni di estremita' nell'asta analizzata
  dNTM([1,3,5],ne)=-dNTM([1,3,5],ne);
end

%%%%%%%%%%%%%%%%%%%%%%  FINE ELABORAZIONE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%  RAPPRESENTAZIONE  GRAFICA  %%%%%%%%%%%%%%%%%%%%%%%%

% Definizione della finestra
figure(1)
clf
hold on

dXmin=min(dXY(:,1));
dXmax=max(dXY(:,1));
dYmin=min(dXY(:,2));
dYmax=max(dXY(:,2));

maxdeltaX=dXmax-dXmin;
maxdeltaY=dYmax-dYmin;

if (maxdeltaY==0)
  maxdeltaY=maxdeltaX;
elseif (maxdeltaX==0)
  maxdeltaX=maxdeltaY;
end

dmax=max([maxdeltaX, maxdeltaY]);  % dimensione massima (X o Y) della struttura

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

% Configurazione iniziale
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);

  plot([dXYA(1),dXYB(1)],[dXYA(2),dXYB(2)],'b-')
  text((dXYA(1)+dXYB(1))/2,(dXYA(2)+dXYB(2))/2-(dLn/20),sprintf('%d',n),'color',[0,0,1]);
end


for nn=1:nNodi
  dXYnn=dXY(nn,:);
  plot([dXYnn(1)-dmax/100,dXYnn(1)+dmax/100],[dXYnn(2),dXYnn(2)],'k-')      
  plot([dXYnn(1),dXYnn(1)],[dXYnn(2)-dmax/100,dXYnn(2)+dmax/100],'k-')      
  text(dXYnn(1)+dmax/30,dXYnn(2)+dmax/30,sprintf('%d',nn))
end



% Configurazione deformata
dUmax=0;
ds=[0:.01:1];
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  dRt=[ dcosAlpha, dsinAlpha;
       -dsinAlpha, dcosAlpha];

  dEAn=dEAI(n,1);
  dEIn=dEAI(n,2);

  dUnXY=dU(nIncGdl(n,:),1);
  dUn=dUnXY;
  dUn([1,2],1)=dRt*dUnXY([1,2],1);
  dUn([4,5],1)=dRt*dUnXY([4,5],1);

  dq=dpqm(n,2);
  dUmax=max([dUmax,abs([(1-3*ds.^2+2*ds.^3)*dUn(2,1)+dLn*(ds-2*ds.^2+ds.^3)*dUn(3,1)+(3*ds.^2-2*ds.^3)*dUn(5,1)+dLn*(-ds.^2+ds.^3)*dUn(6,1)]+dq/dEIn/24*dLn^4*((1-ds).^2).*(ds.^2))]);
end

dUmax

if (dUmax>0)
  dAmplificazione=(dmax/10)/dUmax;  % Lo spostamento massimo verra' rappresentato come dmax/10
else
  dAmplificazione=1;
end    

ds=[0:.01:1];
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  dRt=[ dcosAlpha, dsinAlpha;
       -dsinAlpha, dcosAlpha];

  dEAn=dEAI(n,1);
  dEIn=dEAI(n,2);

  dUnXY=dU(nIncGdl(n,:),1);
  dUn=dUnXY;
  dUn([1,2],1)=dRt*dUnXY([1,2],1);
  dUn([4,5],1)=dRt*dUnXY([4,5],1);

  dXYd=[(1-ds)*dXYA(1)+ds*dXYB(1);
        (1-ds)*dXYA(2)+ds*dXYB(2)];

  dq=dpqm(n,2);
  dXYd=dXYd+dAmplificazione*dRt'*[(1-ds)*dUn(1,1)+ds*dUn(4,1)+(dp/dEAn/2*dLn^2)*ds.*(1-ds);
                                  [(1-3*ds.^2+2*ds.^3)*dUn(2,1)+dLn*(ds-2*ds.^2+ds.^3)*dUn(3,1)+(3*ds.^2-2*ds.^3)*dUn(5,1)+dLn*(-ds.^2+ds.^3)*dUn(6,1)]+dq/dEIn/24*dLn^4*((1-ds).^2).*(ds.^2)];
  plot(dXYd(1,:),dXYd(2,:),'r-')
end



% Rappresentazione azioni interne (azioni assiali) sulla configurazione iniziale 
% Definizione della finestra
figure(2)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

dNmax=max(abs([dNTM(1,:),dNTM(4,:)]));
if (dNmax>0)
  dNAmplificazione=(dmax/10)/dNmax;  % Il diagramma verra' rappresentato con ordinate massime pari a dmax/10
else
  dNAmplificazione=1;
end    

title('Azione assiale')
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  plot([dXYA(1),dXYB(1)],[dXYA(2),dXYB(2)],'k-')
  
  fill([dXYB(1), dXYA(1), dXYA(1)+dNAmplificazione*dNTM(1,n)*dsinAlpha, dXYB(1)+dNAmplificazione*dNTM(4,n)*dsinAlpha, dXYB(1)],...
       [dXYB(2), dXYA(2), dXYA(2)-dNAmplificazione*dNTM(1,n)*dcosAlpha, dXYB(2)-dNAmplificazione*dNTM(4,n)*dcosAlpha, dXYB(2)],'r-')
end
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  if (dNTM(1,n)>=0)
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(1,n)))
  else
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(1,n)))
  end
  if (dNTM(4,n)>=0)
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(4,n)))
  else
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(4,n)))
  end
end




figure(3)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

dTmax=max(abs([dNTM(2,:),dNTM(5,:)]));
if (dTmax>0)
  dTAmplificazione=(dmax/10)/dTmax;  % Il diagramma verra' rappresentato con ordinate massime pari a dmax/10
else
  dTAmplificazione=1;
end    

title('Azione di taglio')
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  plot([dXYA(1),dXYB(1)],[dXYA(2),dXYB(2)],'k-')
  
  fill([dXYB(1), dXYA(1), dXYA(1)+dTAmplificazione*dNTM(2,n)*dsinAlpha, dXYB(1)+dTAmplificazione*dNTM(5,n)*dsinAlpha, dXYB(1)],...
       [dXYB(2), dXYA(2), dXYA(2)-dTAmplificazione*dNTM(2,n)*dcosAlpha, dXYB(2)-dTAmplificazione*dNTM(5,n)*dcosAlpha, dXYB(2)],'r-')
end
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  if (dNTM(2,n)>=0)
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(2,n)))
  else
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(2,n)))
  end
  if (dNTM(5,n)>=0)
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(5,n)))
  else
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(5,n)))
  end
end




figure(4)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

dMmax=0;
for n=1:nAste
  if (dNTM(2,n)*dNTM(5,n)<0)
    dq=dpqm(n,2);
    dMmax=max([dMmax,abs([dNTM(3,n),dNTM(6,n),dNTM(3,n)-dNTM(2,n)^2/2/dq])]);  
  else
    dMmax=max([dMmax,abs([dNTM(3,n),dNTM(6,n)])]);  
  end
end

if (dMmax>0)
  dMAmplificazione=(dmax/10)/dMmax;  % Il diagramma verra' rappresentato con ordinate massime pari a dmax/10
else
  dMAmplificazione=1;
end    

title('Momento flettente')
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  plot([dXYA(1),dXYB(1)],[dXYA(2),dXYB(2)],'k-')

  dq=dpqm(n,2);
  fill([dXYB(1), dXYA(1), dXYA(1)+[0:.01:1]*(dXYB(1)-dXYA(1))+dMAmplificazione*[dNTM(3,n)+dNTM(2,n)*[0:.01:1]*dLn+dq/2*([0:.01:1]*dLn).^2]*dsinAlpha, dXYB(1)],...
       [dXYB(2), dXYA(2), dXYA(2)+[0:.01:1]*(dXYB(2)-dXYA(2))-dMAmplificazione*[dNTM(3,n)+dNTM(2,n)*[0:.01:1]*dLn+dq/2*([0:.01:1]*dLn).^2]*dcosAlpha, dXYB(2)],'r-')
end
for n=1:nAste
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);

  dLn=norm(dXYB-dXYA);
  
  dcosAlpha=(dXYB(1)-dXYA(1))/dLn;
  dsinAlpha=(dXYB(2)-dXYA(2))/dLn;

  if (dNTM(3,n)>=0)
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(3,n)))
  else
    text(dXYA(1)+(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYA(2)+(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(3,n)))
  end
  if (dNTM(6,n)>=0)
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1-(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1+(dLn/20)*dcosAlpha,sprintf('%g',dNTM(6,n)))
  else
    text(dXYB(1)-(dXYB(1)-dXYA(1))*.1+(dLn/20)*dsinAlpha,dXYB(2)-(dXYB(2)-dXYA(2))*.1-(dLn/20)*dcosAlpha,sprintf('%g',dNTM(6,n)))
  end
end





