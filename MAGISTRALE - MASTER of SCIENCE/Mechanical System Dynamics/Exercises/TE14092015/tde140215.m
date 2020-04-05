clear all
close all

m1 = 10;
EJ1 = 5.0e6;
EA1 = 2.6e8;
m2 = 50;
EJ2 = 1.0e8;
EA2 = 1.3e9;
fmax = 20;
coef = 1.5;

om = coef*2*pi*fmax;
Lmax1 = sqrt(pi^2/om*sqrt(EJ1/m1));
Lmax2 = sqrt(pi^2/om*sqrt(EJ2/m2));


[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam]=loadstructure; %ndof or ngdl

% Plot undeformed structure
figure
dis_stru(posiz,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')

% Assembling of elements from .inp file 
[M,K]=assem(incidenze,l,m,EA,EJ,gamma,idb);

% assembling of lumped parameter elements to be done manually
h1 = 0.02;
h4 = 0.01;

% partitioning of matrices 
% variable ndof (number of d.o.f.s) to be defined

MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);



%..............................................
% natural frequencies and modes of vibration
[eigenvectors eigenvalues]=eig(MFF\KFF);
freq=sqrt(diag(eigenvalues))/2/pi;
[dummy,ind]=sort(freq);

AA = [1/(2*2*pi*dummy(1)), dummy(1)*pi;1/(2*2*pi*dummy(4)), dummy(4)*pi];
XX = inv(AA)*[h1;h4];

% alfa and beta values for structural damping to be assigned manually
alfa=XX(1);
beta=XX(2);
C=alfa*M+beta*K;
CFF=C(1:ndof,1:ndof);

% number of modes to be plotted to be assigned manually
nmodes=3;
modes=ind(1:nmodes);

%scaling factor to be assigned manually (same for all plots, can be
%adjusted)
scale_factor=3;

for ii=1:nmodes
    figure
    diseg2(eigenvectors(:,modes(ii)),scale_factor,incidenze,l,gamma,posiz,idb,xy);
    title(['f_',num2str(ii),'=',num2str(freq(modes(ii))),' Hz'])
    if ii<4
        phi(:,ii)=eigenvectors(:,modes(ii));
    end
end

modalMFF = phi'*MFF*phi;
modalKFF = phi'*KFF*phi;
modalCFF = phi'*CFF*phi;


%..............................................
% FRF

% vector of frequency to be defined
vett_f=0:0.01:fmax;
omega=vett_f*2*pi;

% Forcing vector 
F0=zeros(ndof,1);
ndof13 = idb(13,1);
F0(ndof13,1) = 1;
modalF0 = phi'*F0;
for k=1:length(vett_f)
    A=-omega(k)^2*MFF+sqrt(-1)*omega(k)*CFF+KFF;
    xx(:,k)=(A\F0)*sqrt(-1)*omega(k);
    xxx(:,k)=(A\F0)*(-(omega(k))^2);
    
    modalA=-omega(k)^2*modalMFF+sqrt(-1)*omega(k)*modalCFF+modalKFF;
    modalxx(:,k)=(modalA\modalF0)*sqrt(-1)*omega(k);
    modalxxx(:,k)=(modalA\modalF0)*(-(omega(k))^2);
end

%outputs to be defined, script prepared for two outputs
xxmodale = phi*modalxx;
xxxmodale = phi*modalxxx;           % nonhoccapitoooo

ndof1 = idb(13,1);
ndof2 = idb(16,2);

FRF_YppP_F = xxx(ndof2,:);
FRF_YppP_3f = xxxmodale(ndof2,:);
FRF_XpB_F = xx(ndof1,:);
FRF_XpB_3f = xxmodale(ndof1,:);

figure
subplot(211)
plot(vett_f,abs(FRF_XpB_F),'b',vett_f,abs(FRF_XpB_3f),'r')
grid on
title('FRF: Horizontal Velocity of B vs Force')
ylabel('$|X_{pB}/F_B|  [m/(sN)]$','interpreter','latex')
subplot(212)
plot(vett_f,angle(FRF_XpB_F),'b',vett_f,angle(FRF_XpB_3f),'r')
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

figure
subplot(211)
plot(vett_f,abs(FRF_YppP_F),'b',vett_f,abs(FRF_YppP_3f),'r')
grid on
title('FRF: Vertical Acceleration of P vs Force')
ylabel('$|Y_{ppP}/F_B|  [m/(s^2N)]$','interpreter','latex')
subplot(212)
plot(vett_f,angle(FRF_YppP_F),'b',vett_f,angle(FRF_YppP_3f),'r')
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

% distributed load

p0 = 200;
a = gamma(27);
p0_hor = -p0*sin(a);
p0_ver = -p0*cos(a);

F_loc = zeros(ndof,1);
F_glo = zeros(ndof,1);

F_loc(idb(16,:),1)=[p0_hor*(l(27)/2); p0_ver*(l(27)/2); p0_ver*(l(27)^2/12)];
lambda1 = [cos(gamma(27)) sin(gamma(27)) 0; -sin(gamma(27)) cos(gamma(27)) 0; 0 0 1];
F_glo(idb(16,:),1) = F_loc(idb(16,:),1)'*lambda1;

F_loc(idb(15,:),1)=[p0_hor*(l(27)/2+l(28)/2); p0_ver*(l(27)/2+l(28)/2); p0_ver*(-l(27)^2/12+l(28)^2/12)];
lambda2 = [cos(gamma(28)) sin(gamma(28)) 0; -sin(gamma(28)) cos(gamma(28)) 0; 0 0 1];
F_glo(idb(15,:),1) = F_loc(idb(15,:),1)'*lambda2;

F_loc(idb(14,:),1)=[p0_hor*(l(28)/2); p0_ver*(l(28)/2); -p0_ver*(l(28)^2/12)];
F_glo(idb(14,:),1) = F_loc(idb(14,:),1)'*lambda2;

idof_P_ver = idb(16,2);
idof_P_hor = idb(16,1);

xx_stat = inv(KFF)*F_glo;
disp('Vertical Displacement due to distributed load')
X_P_ver_stat = xx_stat(idof_P_ver)
disp('Horizontal Displacement due to distributed load')
X_P_hor_stat = xx_stat(idof_P_hor)

%Uncomment the following lines if you need to analyse the effect of a periodic force
%%..............................................
%% Fourier decomposition of a periodic signal
%% period of signal to be assigned manually
%T=;
%
%% values of signal at equally spaced discrete time values to be assigned
%% manually
%vect_F=[];
%
%fftout=fft(vect_F);
%N=length(vect_F);
%df=1/T;
%fmax=(N/2-1)*df;
%vett_freq=0:df:fmax;
%modf(1)=1/N*abs(fftout(1));
%modf(2:N/2)=2/N*abs(fftout(2:N/2));
%fasf(1:N/2)=angle(fftout(1:N/2));
%
%figure
%subplot 211;bar(vett_freq,modf);grid;
%subplot 212;bar(vett_freq,fasf);grid



