clear all
close all

%% TE 04 09 2016

%% Data
mp = 4;  %[kg/m]
mr = 1;  %[kg/m]
EAp = 4e7;%[N]
EAr = 4e3;%[N]
EJp = 1e5;%[Nm^2]
EJr = 1e3;%[Nm^2]

% control on finite elements length
sf = 2;
fmax = 20; %{Hz]
ommax = 2*pi*fmax*2;

lmaxp = sqrt(pi^2/ommax * sqrt(EJp/mp))
lmaxr = sqrt(pi^2/ommax * sqrt(EJr/mr))

[file_i,xy,nnod,sizee,idb,ndof,incidence,l,gamma,m,EA,EJ,posiz,nbeam]=MeccFEM2_loadstructure('IMP04092016')

% Assembling of elements from .inp file 
[M,K]=assem(incidence,l,m,EA,EJ,gamma,idb);

% assembling of lumped parameter elements to be done manually

% alfa and beta values for structural damping to be assigned manually
alfa=;
beta=;
C=alfa*M+beta*K;

% partitioning of matrices 
% variable ndof (number of d.o.f.s) to be defined
ndof=;
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);
CFF=C(1:ndof,1:ndof);


%..............................................
% natural frequencies and modes of vibration
[eigenvectors eigenvalues]=eig(MFF\KFF);
freq=sqrt(diag(eigenvalues))/2/pi
[dummy,ind]=sort(freq);

% number of modes to be plotted to be assigned manually
nmodes=;
modes=ind(1:nmodes);

%scaling factor to be assigned manually (same for all plots, can be
%adjusted)
scale_factor=;

for ii=1:nmodes
    figure
    diseg2(eigenvalues(:,modes(ii)),scale_factor,incidence,l,gamma,posiz,idb,xy);
    title(['f_',num2str(ii),'=',num2str(fre(modes(ii))),' Hz'])
end

%..............................................
% FRF

% vector of frequency to be defined
vett_f=;
omega=vett_f*2*pi;

% Forcing vector 
F0=;
for k=1:length(vett_f)
    A=-omega(k)^2*MFF+i*omega(k)*CFF+KFF;
    x=A\F0;
    %outputs to be defined, script prepared for two outputs 
    yg1=;
    yg2=;
    mod1(k)=abs(yg1);
    fas1(k)=angle(yg1);
    mod2(k)=abs(yg2);
    fas2(k)=angle(yg2);
end

figure
subplot 211;plot(vett_f,mod1);grid;
subplot 212;plot(vett_f,fas1);grid

figure
subplot 211;plot(vett_f,mod2);grid;
subplot 212;plot(vett_f,fas2);grid


%..............................................
% Fourier decomposition of a periodic signal
% period of signal to be assigned manually
T=;

% values of signal at equally spaced discrete time values to be assigned
% manually
vect_F=[];

fftout=fft(vect_F);
N=length(vect_F);
df=1/T;
fmax=(N/2-1)*df;
vett_freq=0:df:fmax;
modf(1)=1/N*abs(fftout(1));
modf(2:N/2)=2/N*abs(fftout(2:N/2));
fasf(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211;bar(vett_freq,modf);grid;
subplot 212;bar(vett_freq,fasf);grid



