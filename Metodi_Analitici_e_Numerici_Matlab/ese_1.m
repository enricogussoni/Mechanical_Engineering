% Esercizio sui metodi di quadratura

% Punto 1

clear all
close all
clc
figure(1);
a = 0; b = 2*pi;
x_dis = linspace(a,b,1000);
f = @(x) x.*sin(x)./(2*pi);
y_dis=f(x_dis);
plot(x_dis,y_dis, 'r', 'linewidth',3 )
% opzioni grafico
grid on
title('f(x)=xsin(x)')
xlabel('x')
ylabel('f(x)')
%saveas(gcf, 'graficof','epsc2')

%% Punto 4

N = [1:20];
V_ex = -1;
V_PMC = zeros(1,N(end));
V_TC = zeros(1,N(end));
V_SC = zeros(1,N(end));
for i = N
    V_PMC(i) = pmedcomp(a,b,i,f);
    V_TC(i) = trapcomp(a,b,i,f);
    V_SC(i) = simpcomp(a,b,i,f);
end
figure(2);
plot(N,V_PMC,'*',N,V_TC,'o',N,V_SC,'d',N, V_ex * ones(1,20), 'linewidth', 3)
grid on
legend('Pto medio composito','Trapezio composito','Simpson composito','V', 1) 
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Numero sottointervalli')
ylabel('Integrale approssimato')
title('Andamento dell` integrale approssimato')
%saveas(gcf, 'ConfrontoIntegrali','epsc2')

%% Punto 5

err_PMC= abs(V_ex - V_PMC);
err_TC= abs(V_ex - V_TC);
err_SC=abs(V_ex - V_SC);
H = (b-a)./N;
figure(3);
loglog(H,err_PMC,'*',H,err_TC,'o',H,err_SC,'d',H,H.^2,H,H.^4, 'linewidth', 3)
legend('err-PMC','err-TC','err-SC','H^2','H^4', 2)
% opzioni grafico
set(gca,'FontSize',20)
xlabel('Ampiezza dei sottointervalli')
ylabel('Errore')
title('Errore di quadratura in scala logaritmica')
%saveas(gcf, 'ConfrontoErrori','epsc2')

rapp = err_PMC ./err_TC;
disp(rapp(end-4:end))

%% Punto 6

toll=1e-5;
d2f_int = @(x) (2*cos(x)-x.*sin(x))./(2*pi);
K2 = max(abs(d2f_int(x_dis)));
Npmc = ceil(sqrt( ((b-a)^3 *K2)/(24*toll) ) )

VPM = pmedcomp(a,b,Npmc,f);
errPMCvero = abs(V_ex - VPM)

