function DL01 = precarico1()

global G 
global m1 m2 m3 k1
global theta0 ddy1 ddy2 ddy3 DL ddDL 

[dy1,ddy1,dy2,ddy2,dy3,ddy3,DL,dDL,ddDL] = energia_potenziale1(theta0);

% Allungamento statico dovuto al precarico
DL01 = -(G*(m1*dy1+m2*dy2+m3*dy3)/(k1*dDL)); 