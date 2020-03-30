% Lab 5 Es 4
%INCOMPLETO

sigma   = [0.18   0.3   0.5    0.6    0.72  0.75   0.8   0.9   1.0   ];
epsilon = [0.0005 0.001 0.0013 0.0015 0.002 0.0045 0.006 0.007 0.0085];

% figure(1)
% plot(sforzo,def)
% xlabel('sforzo in 1000 Kgf/cm^2')
% ylabel('deformazione in cm/cm')

%1
x_sigma=linspace(sigma(1),sigma(end),1000);

% Polinomio di Lagrange
grado=length(sigma)-1;
p=polyfit(sigma,epsilon,grado);
def_p=polyval(p,x_sigma);

figure(2)
plot(x_sigma,def)
xlabel('sforzo')
ylabel('deformazione')
title('Curva sforzo-deformazione')

% Interpolazione composita lineare
def_l=interp1(sigma,epsilon,x_sigma);

figure(3)
plot(x_sigma,def_l)
xlabel('sforzo')
ylabel('deformazione')
title('Grafico sforzo-deformazione')

%2
sigma_sp=[0.4 0.65];
epsilon_p=polyval(def_p,sigma_sp);
epsilon_l=interp1(sigma,epsilon,sigma_sp);