function [clcd]=Re_medio(clRem,clReM,cdRem,cdReM,Re,ReM,Rem)

% Interpolazione lineare di CL e CD a Re conoscendo i valori per Re_massimo (ReM)
% e Re_minimo (Rem)

clRe=((clReM-clRem)./(ReM-Rem)).*(Re-Rem)+clRem;
cdRe=((cdReM-cdRem)./(ReM-Rem)).*(Re-Rem)+cdRem;
clcd=[clRe;cdRe];
end
