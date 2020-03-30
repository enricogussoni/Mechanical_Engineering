% Lastra piana perpendicolare al flusso

L=0.8;                  % larghezza piastra frenante
b=0.4;                  % altezza " "
A=L*b;                  % superficie frontale " "

v=100/3.6;              % velocità inizio frenata
patm=1e5;               % pressione atmosferica
pm=b/2*9.81*1000;       % pressione in acqua
p=v^2*1000/2;           % pressione spingente

F=p*A;                  % forza frenante

