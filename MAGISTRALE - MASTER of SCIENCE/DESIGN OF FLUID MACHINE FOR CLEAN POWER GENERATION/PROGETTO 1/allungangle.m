function [CL, CD, angle] = allungangle(CL,CD,angle)

% Estende artificialmente le curve di CL e CD oltre il limite negativo dei
% dati. Approssimazione per il Pitch-to-Feather.


pos0 = find(abs(CL-0)<1e-2);
pos0 = pos0(1);
aaa = angle(pos0);
postall = find(angle==12);
CL=[-CL(postall:-1:pos0),CL(pos0+1:end)];
CD=[CD(postall:-1:pos0),CD(pos0+1:end)];
angle = (2*aaa-angle(postall)):0.1:angle(end);