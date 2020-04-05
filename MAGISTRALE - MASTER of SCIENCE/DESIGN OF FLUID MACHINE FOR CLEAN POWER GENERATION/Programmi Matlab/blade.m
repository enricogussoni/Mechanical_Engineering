function [matl,matd]=blade(clh,cdh,clp,cdp,clt,cdt,rh,rp,rt,r)
r1=r(r<=rp);
f=length(clh);
onel=zeros(f,r1);
oned=onel;
for i=1:f
    for j=1:length(r1)
        onel(i,j)=((clp(i)-clh(i))/(rp-rh))*(r1(j)-rh)+clh(i);
        oned(i,j)=((cdp(i)-cdh(i))/(rp-rh))*(r1(j)-rh)+cdh(i);
    end
end
r2=r(r>rp);
twol=zeros(f,r2);
twod=twol;
for i=1:f
    for j=1:length(r2)
        twol(i,j)=((clt(i)-clp(i))/(rt-rp))*(r2(j)-rp)+clp(i);
        twod(i,j)=((cdt(i)-cdp(i))/(rt-rp))*(r2(j)-rp)+cdp(i);
    end
end
matl=[onel,twol];
matd=[oned,twod];
end

