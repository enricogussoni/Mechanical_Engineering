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

[file_i,xy,nnod,sizee,idb,ndof,incidence,l,gamma,m,EA,EJ,posiz,nbeam]=MeccFEM2_loadstructure('040915')

% Plot of undeformed structure
MeccFEM2_plotStructure(posiz,l,gamma,xy)

% Assembling of elements from .inp file 
[M,K]=assem(incidence,l,m,EA,EJ,gamma,idb);

% Partitioning
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);

MFC = M(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

%% Natural frequencies and modes of vibration
[modes Om2]=eig(MFF\KFF);
freq=sqrt(diag(Om2))/2/pi
[freqord,ordmode]=sort(freq);
om_dis=sqrt(Om2);
% om=sort(om_dis);

% number of modes to be plotted to be assigned manually
nmodes=3;

%scaling factor to be assigned manually (same for all plots, can be
%adjusted)
scale_factor=2;

scaleFactor=2;        % a scelta -> ampiezza del disegno 
for ii=1:3
    mode=modes(:,ordmode(ii));
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,posiz,idb,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(freqord(ii))])
end  


% alfa and beta values for structural damping to be assigned manually
h1 = 0.02;
h2 = 0.03;
h = [h1,h2];

A = [1/(2*2*pi*freqord(1)) pi*freqord(1) ; 1/(2*2*pi*freqord(2)) pi*freqord(2)];
x = A^(-1)*h';

alfa=x(1);
beta=x(2);
C=alfa*M+beta*K;

RFF = C(1:ndof,1:ndof);

% FRF
freq = linspace(0,fmax,1000);
om = 2*pi*freq;
f = zeros(1,ndof);
ngdlF = idb(10,2);
f(ngdlF)=1;

for k=1:length(freq)
    xx(:,k)=inv(-om(k)^2*MFF+sqrt(-1)*om(k)*RFF+KFF)*f';
end

ngdlj = idb(17,2);
FRFj = xx(ngdlj,:);

figure
subplot(2,1,1)
plot(freq, abs(FRFj))
subplot(2,1,2)
plot(freq, angle(FRFj))


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



