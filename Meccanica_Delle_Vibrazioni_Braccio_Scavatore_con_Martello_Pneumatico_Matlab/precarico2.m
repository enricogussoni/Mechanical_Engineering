function DL20 = precarico2()

global G 
global m1 m2 m3 k2
global theta10 theta20

[dy4,ddy4,dy5,ddy5,dy6,ddy6,DL2,dDL2,ddDL2] = energia_potenziale2(theta10,theta20);

% Allungamento statico dovuto al precarico
DL20 = -(G*(m1*dy4+m2*dy5+m3*dy6)/(k2*dDL2)); 