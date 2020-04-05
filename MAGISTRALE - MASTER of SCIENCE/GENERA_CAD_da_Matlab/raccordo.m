clc
clear all
close all

M=load('naca.txt');
x=M(1:end,1);
y=M(1:end,2);
z=[0 1];
f = [1 2 ];
% costruzione matrici

for j=1:length(z)
    
  X(:,j) = x*f(j); % moltiplicato per j cosi per cambiare dimensione
  Y(:,j) = y*f(j);
  Z(:,j) = z(j)*ones(length(x),1);    
end

[srf] = crea_sup(X,Y,Z,'raccordo.igs');

