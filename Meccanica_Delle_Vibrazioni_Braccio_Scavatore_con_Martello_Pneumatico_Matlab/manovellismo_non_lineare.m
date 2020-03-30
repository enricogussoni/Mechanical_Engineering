function dy = manovellismo_non_lineare(t,y)

global G 
global m1 m2 m3 k1 r1
global z0 wz psiz
global DL01 ddy1 ddy2 ddy3 ddDL

thetap = y(1);
theta = y(2);

z = z0*sin(wz*t+psiz);

M = energia_cinetica();
[dy1,ddy1,dy2,ddy2,dy3,ddy3,DL,dDL,ddDL] = energia_potenziale1(theta);

dV = k1*(DL+DL01)*dDL+m1*G*dy1+m2*G*dy2+m3*G*dy3;
dD = r1*dDL^2*thetap;
Q=z*dy3;

dy = [-(dD+dV-Q)/M; thetap];
