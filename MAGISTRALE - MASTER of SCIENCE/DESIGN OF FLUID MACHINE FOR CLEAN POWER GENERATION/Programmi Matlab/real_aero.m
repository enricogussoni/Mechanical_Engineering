function [mat]=real_aero(alfa,cl,cd,mioalfa,ar)

% Aspect ratio correction for airfoil properties

clmio=interp1(alfa,cl,mioalfa);
cind=clmio./(pi*ar);
alfar=mioalfa-abs(cind)*180/pi;
if sum(alfar<-180)>0
    for i=1:length(alfar)
        if alfar(i)>-180
            appog(i)=alfar(i);
        else
            appog(i)=alfar(i)+360;
        end
    end
else
    appog=alfar;
end
clr=interp1(alfa,cl,appog);
cdr=interp1(alfa,cd,appog)+(clr.^2)./(pi*ar);
mat=[alfar,clr,cdr];
end
