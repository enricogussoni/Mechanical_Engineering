clear all;
close all;

%% Punto1

t=5.8;
w=100;
a=0.05*t;
c=5*a;
phiA=pi/2;
phiC=0;
P=40;
Dout=230;
D=Dout-t;
St=P*D/(2*t);
Sb=0;

Ka = SIF(t,w,a,c,phiA,St,Sb)/sqrt(1000)
Kc = SIF(t,w,a,c,phiC,St,Sb)/sqrt(1000)

KIC=200;
if Ka<KIC && Kc<KIC
    disp('No brittle fracture');
end
%% Punto 2
P=30;
UTS=1000;
md=1;
ma=0.7;
sige=0.4*UTS;
kf=1;

sigA=P*D/(4*t);
sigR=-P;
sigC=P*D/(2*t);

sigAm=sigA/2;  %Alternati uguali, non definisco
sigRm=sigR/2;
sigCm=sigC/2;

sines=sqrt(sigAm^2+sigCm^2+sigRm^2-sigAm*sigCm-sigRm*sigCm-sigRm*sigAm)
In=sigAm+sigRm+sigCm;
sige_prime=md*ma*sige/kf;
eta=sige_prime*(1-In/UTS)/sines

% Finite life assessment

S1=0.9*UTS;
S2=sige_prime;
N1=10^3;
N2=5*10^6;
Nmin=12000;

b=log(S2/S1)/log(N2/N1);
A=S1/N1^b;

N=(sines/(A*(1-In/UTS)))^(1/b)

if N>Nmin
    disp('Minimum number of cycle satisfied');
end

%% Point 3
C = 2.53e-13;
m =3;
St=P*D/(2*t);
a = 0.29;
c = 1.45;
a_f = 0.95*t;
N = 1;
K_c = 6324.55; % unstable crack propagation stress intensity factor
a_v = [a];
c_v = [c];
KA = SIF(t,w,a,c,pi/2,St,0);
KC = SIF(t,w,a,c, 0,St,0);

while a<a_f % enters only if crack lower than 95% thickness
N=N+1;
K_A = SIF(t,w,a,c,pi/2,St,0); % actual critical value for A
K_C = SIF(t,w,a,c, 0,St,0); % actual critical value for C
if (K_A>=K_c)|(K_C>=K_c) % Let's check neither of these 2
break
end

a = a+C*K_A^m;
c = c+C*K_C^m;
a_v = [a_v, a];
c_v = [c_v, c];
KA = [KA, K_A];
KC = [KC, K_C];
end

if (K_A<K_c)&&(K_C<K_c)
disp('95% thickness reached after a number of cycles:');
N
else
disp('unstable propagation reached after a number of cycles:');
N
end
a_final_value = a
c_final_value = c

figure(1)
plot([1:N],a_v,[1:N],c_v);
legend('depth a','semiaxis c')
xlabel('cycles');
ylabel('length [mm]');
title('crack growth');
figure(2)
plot([1:N],KA,[1:N],KC); hold on; plot([1 N],[K_c K_c],'r');
legend('SIF at a','SIF at c','K_{IC}')
xlabel('cycles');
ylabel('magnitude [MPa$$\sqrt{mm}$$]','Interpreter','latex','FontSize',13);
title('SIF growth');

%% Point 4
P=40;
rm=D/2;
a=c;
St=P*D/(2*t);

M=sqrt(1+3.2*a^2/(2*rm*t));
fw=1/(sqrt(cos(pi*a/w)));
a=a*10^(-3);
K=M*fw*St*sqrt(pi*a);

if K<KIC
    disp('Static assessment verified')
end




