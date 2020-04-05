function [matl,matd]=blade(clh,cdh,clp,cdp,clt,cdt,rh,rp,rt,r)

% Interpola linearmente i valori di CL e CD su tutta la pala a partire da quelli su
% HUB, PRIMARY e TIP (già interpolati su Re)
%output: [CL,CD] matrici n*m con n=#alpha, m=#raggi

r1=r(r<=rp);
ff=length(clh);% # alpha (o beta_inf che dir si voglia)
onel=zeros(ff,length(r1));
oned=onel;
for i=1:ff
    for j=1:length(r1)
        onel(i,j)=((clp(i)-clh(i))/(rp-rh))*(r1(j)-rh)+clh(i);
        oned(i,j)=((cdp(i)-cdh(i))/(rp-rh))*(r1(j)-rh)+cdh(i);
    end
end
r2=r(r>rp);
twol=zeros(ff,length(r2));
twod=twol;
for i=1:ff
    for j=1:length(r2)
        twol(i,j)=((clt(i)-clp(i))/(rt-rp))*(r2(j)-rp)+clp(i);
        twod(i,j)=((cdt(i)-cdp(i))/(rt-rp))*(r2(j)-rp)+cdp(i);
    end
end
matl=[onel,twol];
matd=[oned,twod];
end

