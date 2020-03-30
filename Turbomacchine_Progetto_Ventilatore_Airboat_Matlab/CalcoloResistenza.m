%% Accelerazione
% http://progettazionenautica.blogspot.it/2012/03/calcolo-della-resistenza-dattrito.html
v_max=28.28; % m/s
dens=1000; % kg/m^3
lunghezza=5; % m (bagnata)
larghezza=2.2; % m
S=lunghezza*larghezza;
visc_cin=834.2e-6; % m^2/s (tab. viscosità Marco Landonio)

Re=v_max*lunghezza/visc_cin; % Reynolds
cf=0.075/(log(Re)-2)^2;

v_vect=0:0.56:28.28;
resistenza= cf*0.5*S*dens*v_vect.^2;

figure(1)
title('Andamento resistenza')
plot(v_vect,resistenza);
xlabel('Velocità [m/s]')
ylabel('Resistenza [N]')
