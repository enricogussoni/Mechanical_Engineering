

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Corso di Metodi di Calcolo delle Strutture       %
%            Analisi statica di problemi piani           %
%                     25/01/2016                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
clc
close all


bProb=0;  % =0 per problema piano negli sforzi
          % =1 per problema piano nelle deofrmazioni

bNumNEF=0; % =0 per non indicare i numeri degli EF nel disegno della deformata
           % =1 per indicare i numeri degli EF nel disegno della deformata

bMens=1; % =0 per analizzare il problema accademico
         % =1 per analizzare la mensola con carico concentrato all'estremita'


if (bMens==1) %analisi mensola 
  % Definizione delle coordinate dei nodi della struttura
  dL=4000; %[mm] - lunghezza trave
  dh=200; %[mm] - altezza trave

  nn=21;
  nNodiX=nn*10+1;
  nNodiY=nn;

  dXY=zeros([nNodiX*nNodiY,2]);

  for ny=1:nNodiY
    for nx=1:nNodiX
      nn=nx+(ny-1)*nNodiX;
      dXY(nn,1)=(nx-1)*dL/(nNodiX-1);
      dXY(nn,2)=(ny-1)*dh/(nNodiY-1);
    end
  end  
  nNodi=size(dXY,1);

  nGdlTot=2*nNodi;


  % Definizione della tabella delle incidenze dei nodi
  nIncNodi=zeros([(nNodiX-1)*(nNodiY-1)*2,3]);

  ne=0;
  for ny=1:nNodiY-1
    for nx=1:nNodiX-1
      ne=ne+1;
      nIncNodi(ne,:)=[(ny-1)*nNodiX+nx, (ny-1)*nNodiX+nx+1, ny*nNodiX+nx];
      ne=ne+1;
      nIncNodi(ne,:)=[(ny-1)*nNodiX+nx+1, ny*nNodiX+nx+1, ny*nNodiX+nx];
    end
  end  
  nEF=size(nIncNodi,1);

  % Definizione della tabella delle incidenze dei gradi di liberta'
  nIncGdl=[nIncNodi(:,1)*2-1, nIncNodi(:,2)*2-1, nIncNodi(:,3)*2-1, nIncNodi(:,1)*2, nIncNodi(:,2)*2, nIncNodi(:,3)*2];

   
  % Assegnazione dei parametri (E,ni,t)

  %E=[N/mm2]
  %ni=[]
  %t=[mm]

  dEnit=[206000*ones([nEF,1]), 0.3*ones([nEF,1]), 10*ones([nEF,1])];  


 
 
  % Assegnazione dei gradi di liberta' e dei vincoli del sistema
  nUv=sort([[1:2*nNodiX:2*(nNodiX*nNodiY)],[2:2*nNodiX:2*(nNodiX*nNodiY)]]);
  nUl=[1:nGdlTot];
  nUl(nUv)=[];

  % Assegnazione dei carichi
  dF=zeros([nGdlTot,1]);

  %Carichi concentrati
  dF(2*(nNodiX*nNodiY),1)=-5000;  %[N]

  %Carichi distribuiti
  %     [  fxAB,fyAB, fxBC,fyBC, fxCA,fyCA]
  dfsup=zeros([nEF,6]);

  %     [ bx, by]
  dfvol=zeros([nEF,2]);
  
  % Assegnazione degli spostamenti nei nodi vincolati
  dUv=zeros(size(nUv'));

else  %analisi problema accademico
  dL=1000; %[mm]

  %Coordinate dei nodi
  dXY=[ 0, 0;
       dL, 0;
       dL, dL;
        0, dL];

  nNodi=size(dXY,1);  %Numero totale dei nodi
  


  %Elementi finiti (EF) e connessioni (tabelle delle incidenze)
  nIncNodi=[1,2,4;
            2,3,4];

  %nIncGdl=[1, 3, 7, 2, 4, 8;
  %         3, 5, 7, 4, 6, 8]
  nIncGdl=[nIncNodi(:,1)*2-1, nIncNodi(:,2)*2-1, nIncNodi(:,3)*2-1, nIncNodi(:,1)*2, nIncNodi(:,2)*2, nIncNodi(:,3)*2];

  nEF=size(nIncNodi,1);

  nGdlTot=2*nNodi;   %Numero totale dei gradi di liberta' (gdl)



  %Parametri del materiale e spessore (Modulo di Young "E", coefficiente di Poisson "ni", spessore "t")
  %    [    E,    ni,     t]
  dEnit=[206000,   0.3,   10;
         206000,   0.3,   10];

    

  %Gdl da vincolare 
  nUv=[1,2,4,7]; 

  %Gdl da non vincolare
  %nUl=[3,5,6,8]; 
  nUl=[1:nGdlTot];
  nUl(nUv)=[];


  %Vettore termini noti
  dF=zeros([nGdlTot,1]);

  %Carichi distribuiti
  %     [  fxAB,fyAB, fxBC,fyBC, fxCA,fyCA]
  dfsup=[     0,   0,    0,   0,    0,   0;
            200,   0,    0,   0,    0,   0];

  %     [ bx, by]
  dfvol=[  0,  0;
           0,  0];
   
  %Cedimenti vincolari assegnati
  dUv=[0,0,0,0]';
end

%%%%%%%%%%%%%%%%%%%%%%% Termine assegnazione dati %%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Elaborazione %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Costruzione della matrice di rigidezza globale (K)
%dK=zeros([nGdlTot,nGdlTot]);
dK=spalloc(nGdlTot,nGdlTot,36);
for ne=1:nEF
  nA=nIncNodi(ne,1);  %nA=numero del primo nodo dell'EF
  nB=nIncNodi(ne,2);  %nB=numero del secondo nodo dell'EF
  nC=nIncNodi(ne,3);  %nB=numero del terzo nodo dell'EF

  dXYA=dXY(nA,:);  %dXYA=coordinate del primo nodo dell'EF
  dXYB=dXY(nB,:);  %dXYB=coordinate del secondo nodo dell'EF
  dXYC=dXY(nC,:);  %dXYC=coordinate del terzo nodo dell'EF

  dJe=[1,dXYA;
       1,dXYB;
       1,dXYC];
  dAe=abs(det(dJe))/2;
  
  dEne=dEnit(ne,1); %Modulo di Young dell'EF analizzato (numero "ne")
  dnine=dEnit(ne,2); %Coefficiente di Poisson dell'EF analizzato (numero "ne")
  dtne=dEnit(ne,3); %Spessore dell'EF analizzato (numero "ne")

  if (bProb)
    dEne=dEne/(1-dnine^2); %parametro E  per problema piano nelle deformazioni
    dnine=dnine/(1-dnine); %parametro ni per problema piano nelle deformazioni
  end
  dEmat=dEne/(1-dnine^2)*[    1, dnine,           0;
                          dnine,     1,           0;
                              0,     0, (1-dnine)/2];

  dBe=[dXYB(2)-dXYC(2), dXYC(2)-dXYA(2), dXYA(2)-dXYB(2),               0,               0,               0;
                     0,               0,               0, dXYC(1)-dXYB(1), dXYA(1)-dXYC(1), dXYB(1)-dXYA(1);
       dXYC(1)-dXYB(1), dXYA(1)-dXYC(1), dXYB(1)-dXYA(1), dXYB(2)-dXYC(2), dXYC(2)-dXYA(2), dXYA(2)-dXYB(2)]/det(dJe);
  
 
  %Matrice di rigidezza della trave analizzata
  dk=dtne*dAe*dBe'*dEmat*dBe; 
  
  df=zeros([6,1]);
  dLca=norm(dXYC-dXYA);  %dL=lunghezza del tratto CA
  dLab=norm(dXYA-dXYB);  %dL=lunghezza del tratto AB
  dLbc=norm(dXYB-dXYC);  %dL=lunghezza del tratto BC

% dfsup(ne,:)=[ fxAB,fyAB, fxBC,fyBC, fxCA,fyCA]
  df=dtne/2*[dLca*dfsup(ne,5)+dLab*dfsup(ne,1);
             dLab*dfsup(ne,1)+dLbc*dfsup(ne,3);
             dLbc*dfsup(ne,3)+dLca*dfsup(ne,5);
             dLca*dfsup(ne,6)+dLab*dfsup(ne,2);
             dLab*dfsup(ne,2)+dLbc*dfsup(ne,4);
             dLbc*dfsup(ne,4)+dLca*dfsup(ne,6)];
  
% dfvol(ne,:)=[ bx, by]
  df=df+dAe*dtne/3*[dfvol(ne,1);
                    dfvol(ne,1);
                    dfvol(ne,1);
                    dfvol(ne,2);
                    dfvol(ne,2);
                    dfvol(ne,2)];
  
  
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




%%%%%%%%%%%%%%% Determinazione degli sforzi %%%%%%%%%%%%%%%%%%%
%dSigma(ne,:)=[sigma_xx, sigma_yy, tau_xy]
dSigma=zeros([nEF,3]);
for ne=1:nEF
  nA=nIncNodi(ne,1);  %nA=numero del primo nodo dell'EF
  nB=nIncNodi(ne,2);  %nB=numero del secondo nodo dell'EF
  nC=nIncNodi(ne,3);  %nB=numero del terzo nodo dell'EF

  dXYA=dXY(nA,:);  %dXYA=coordinate del primo nodo dell'EF
  dXYB=dXY(nB,:);  %dXYB=coordinate del secondo nodo dell'EF
  dXYC=dXY(nC,:);  %dXYC=coordinate del terzo nodo dell'EF

  dJe=[1,dXYA;
       1,dXYB;
       1,dXYC];
  
  
  dEne=dEnit(ne,1); %Modulo di Young dell'EF analizzato (numero "ne")
  dnine=dEnit(ne,2); %Coefficiente di Poisson dell'EF analizzato (numero "ne")

  if (bProb)
    dEne=dEne/(1-dnine^2); %parametro E  per problema piano nelle deformazioni
    dnine=dnine/(1-dnine); %parametro ni per problema piano nelle deformazioni
  end
  dEmat=dEne/(1-dnine^2)*[    1, dnine,           0;
                          dnine,     1,           0;
                              0,     0, (1-dnine)/2];

  dBe=[dXYB(2)-dXYC(2), dXYC(2)-dXYA(2), dXYA(2)-dXYB(2),               0,               0,               0;
                     0,               0,               0, dXYC(1)-dXYB(1), dXYA(1)-dXYC(1), dXYB(1)-dXYA(1);
       dXYC(1)-dXYB(1), dXYA(1)-dXYC(1), dXYB(1)-dXYA(1), dXYB(2)-dXYC(2), dXYC(2)-dXYA(2), dXYA(2)-dXYB(2)]/det(dJe);
  
 
  %Assemblaggio della matrice dk e del vettore df
  nv=nIncGdl(ne,:); %Vettore dei gdl relativi ai nodi dell'asta analizzata
  dUabc=dU(nv,1);
  ds=dEmat*dBe*dUabc;
  dSigma(ne,:)=ds';
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

dmax=max([maxdeltaX, maxdeltaY]);  % dimensione massima (X o Y) della struttura

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

% Configurazione iniziale
for n=1:nEF
  nA=nIncNodi(n,1);
  nB=nIncNodi(n,2);
  nC=nIncNodi(n,3);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);
  dXYC=dXY(nC,:);

  plot([dXYA(1),dXYB(1),dXYC(1),dXYA(1)],[dXYA(2),dXYB(2),dXYC(2),dXYA(2)],'b-')
  if (bNumNEF)
    text((dXYA(1)+dXYB(1)+dXYC(1))/3,(dXYA(2)+dXYB(2)+dXYC(2))/3,sprintf('%d',n),'color',[0,0,1]);
  end
end





% Configurazione deformata

dUmax=max(abs(dU)); % Spostamento massimo (in direzione X o Y)

if (dUmax>0)
  dAmplificazione=(dmax/10)/dUmax;  % Lo spostamento massimo verra' rappresentato
                                    %    come dmax/10
else
  dAmplificazione=1;
end    
    
for ne=1:nEF
  nA=nIncNodi(ne,1);
  nB=nIncNodi(ne,2);
  nC=nIncNodi(ne,3);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);
  dXYC=dXY(nC,:);

  dXYAd=dXYA+dU(nIncGdl(ne,[1,4]))'*dAmplificazione;
  dXYBd=dXYB+dU(nIncGdl(ne,[2,5]))'*dAmplificazione;
  dXYCd=dXYC+dU(nIncGdl(ne,[3,6]))'*dAmplificazione;
  
  plot([dXYAd(1),dXYBd(1),dXYCd(1),dXYAd(1)],[dXYAd(2),dXYBd(2),dXYCd(2),dXYAd(2)],'r-')      
end




% Rappresentazione sforzi sulla configurazione iniziale 

% Ricerca del massimo sforzo
dsMax=max(max(dSigma));
dsMin=min(min(dSigma));

% Definizione della finestra
figure(2)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

title('s_x_x')
for ne=1:nEF
  nA=nIncNodi(ne,1);
  nB=nIncNodi(ne,2);
  nC=nIncNodi(ne,3);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);
  dXYC=dXY(nC,:);

  fill([dXYA(1), dXYB(1), dXYC(1), dXYA(1)],...
       [dXYA(2), dXYB(2), dXYC(2), dXYA(2)],dSigma(ne,1))
end

caxis([dsMin,dsMax]);
RevJet=jet;
RevJet=RevJet([size(jet,1):-1:1],:);
colormap(RevJet)
colorbar





% Definizione della finestra
figure(3)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

title('s_y_y')
for ne=1:nEF
  nA=nIncNodi(ne,1);
  nB=nIncNodi(ne,2);
  nC=nIncNodi(ne,3);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);
  dXYC=dXY(nC,:);

  fill([dXYA(1), dXYB(1), dXYC(1), dXYA(1)],...
       [dXYA(2), dXYB(2), dXYC(2), dXYA(2)],dSigma(ne,2))
end

caxis([dsMin,dsMax]);
colormap(RevJet)
colorbar





% Definizione della finestra
figure(4)
clf
hold on

axis([dXmin-maxdeltaX/5, dXmax+maxdeltaX/5, dYmin-maxdeltaY/5, dYmax+maxdeltaY/5]);
axis equal

title('t_x_y')
for ne=1:nEF
  nA=nIncNodi(ne,1);
  nB=nIncNodi(ne,2);
  nC=nIncNodi(ne,3);

  dXYA=dXY(nA,:);
  dXYB=dXY(nB,:);
  dXYC=dXY(nC,:);

  fill([dXYA(1), dXYB(1), dXYC(1), dXYA(1)],...
       [dXYA(2), dXYB(2), dXYC(2), dXYA(2)],dSigma(ne,3))
end

caxis([dsMin,dsMax]);
colormap(RevJet)
colorbar





