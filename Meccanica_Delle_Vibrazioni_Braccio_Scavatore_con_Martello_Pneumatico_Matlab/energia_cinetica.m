
function M = energia_cinetica()

global a s theta0 p x
global m1 J1 m2 J2 m3 J3

b2=sqrt(((a+s)*cos(theta0)+p*cos(theta0-(pi*10)/24))^2+((a+s)*sin(theta0)+p*sin(theta0-(pi*10)/24))^2);
b3=sqrt(((a+s)*cos(theta0)+(p+x)*cos(theta0-(pi*10)/24))^2+((a+s)*sin(theta0)+(p+x)*sin(theta0-(pi*10)/24))^2);

M = m1*(a^2)+J1+ m2*(b2^2)+J2+ m3*(b3^2)+J3;

