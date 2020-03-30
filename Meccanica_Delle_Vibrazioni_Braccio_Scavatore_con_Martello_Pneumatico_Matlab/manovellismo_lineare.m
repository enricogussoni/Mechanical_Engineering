function dy = manovellismo_lineare(t,y)

global G 
global m1 m2 m3 k1 r1
global z0 wz psiz
global theta0 DL01 dy1 dy2 DL

thetap = y(1); 
theta = y(2);

z = z0*sin(wz*t+psiz) 

M = energia_cinetica();
[dy1,ddy1,dy2,ddy2,dy3,ddy3,DL,dDL,ddDL] = energia_potenziale1(theta0);
dV = (k1*dDL^2+k1*DL01*ddDL+m1*G*ddy1+m2*G*ddy2+m3*G*ddy3)*(theta-theta0);
dD = r1*dDL^2*thetap;
Q=z*dy3;

dy = [-(dD+dV-Q)/M;
      thetap];
