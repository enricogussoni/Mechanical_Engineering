function [mat]=vit_corr(astall,clstall,cdstall,ar,alfa)
a=astall*pi/180;
al=alfa*pi/180;
if ar<=50
    cdmax=1.11+0.018*ar;
else
    cdmax=2.01;
end
Kl=(clstall-cdmax*sin(a)*cos(a))*sin(a)/((cos(a))^2);
Kd=(cdstall-cdmax*((sin(a))^2))/cos(a);
cl=cdmax*sin(al*2)/2+Kl*((cos(al)).^2)./sin(al);
cd=cdmax*((sin(al)).^2)+Kd*cos(al);
mat=[cl',cd'];
end
