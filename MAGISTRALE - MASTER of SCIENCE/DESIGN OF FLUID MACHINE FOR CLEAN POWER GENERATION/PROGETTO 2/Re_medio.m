function [clcd]=Re_medio(clRem,clReM,cdRem,cdReM,Re,ReM,Rem)
clRe=((clReM-clRem)./(ReM-Rem)).*(Re-Rem)+clRem;
cdRe=((cdReM-cdRem)./(ReM-Rem)).*(Re-Rem)+cdRem;
clcd=[clRe;cdRe];
end
