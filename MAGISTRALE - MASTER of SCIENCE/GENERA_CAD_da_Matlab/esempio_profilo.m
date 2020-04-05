clc
clear all
close all


M=load('naca.txt');
x=M(1:end,1);
y=M(1:end,2);
% z=[0 0.1 0.5 0.8 1];
z=[0 0.1 0.5 0.8 1];
f = [1e-3 2 3 2 1];
% costruzione matrici

% ciclo sulle 3 sezioni a z=[0 10 100] con z coordinate lungo l'aperura

for j=1:length(z)
    
  X(:,j) = x*f(j); % moltiplicato per j cosi per cambiare dimensione
  Y(:,j) = y*f(j);
  Z(:,j) = z(j)*ones(length(x),1);    
end


[srf] = crea_sup(X,Y,Z,'test.igs');

