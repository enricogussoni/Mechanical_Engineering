clear
close all
clc

%% Structure properties
Amon = 52.79e-4;
Amen = 17.05e-4;
Jmon = 4696e-8;
Jmen = 392.9e-8;
E = 2.06e11;
EAmon = E*Amon;
EAmen = E*Amen;
EJmon = E*Jmon;
EJmen = E*Jmen;
rho = 7850;
mmon = rho*Amon;
mmen = rho*Amen;


%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE20072011');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);
ndoc = size(M,1)-ndof;

% Additional mass
Mc = 5;
idfMc = idf(8,[1 2]);
M(idfMc,idfMc)=M(idfMc,idfMc)+[Mc 0; 0 Mc];

% Additional springs
ka = 6.8e7;
Kaloc = [ka 0 0 -ka 0 0;
         0  0 0  0  0 0;
         0  0 0  0  0 0;
         -ka  0 0 ka  0 0;
         0  0 0  0  0 0;
         0  0 0  0  0 0];
gammaa = 2.55;
lambda = [cos(gammaa) sin(gammaa) 0;
         -sin(gammaa) cos(gammaa) 0;
          0           0           1];
Lambda = [lambda zeros(3,3);
          zeros(3,3) lambda];
Kaglob = Lambda'*Kaloc*Lambda;
idfka = [idf(5,:) idf(7,:)];
K(idfka,idfka) = K(idfka,idfka) +Kaglob;

kc = 1e5;
Kc = [kc -kc;-kc kc];
idfkc = [idf(7,2) idf(8,2)];
K(idfkc,idfkc)=K(idfkc,idfkc)+Kc;
      
%% Partitioning of "FF" Mass and Stiffness Matrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes, Om2] = eig(MFF\KFF);
om2 = diag(Om2);
frq = sqrt(om2)/2/pi;

% Sort in ascending order frequencies and mode shapes
[frqord,ordmode] = sort(frq);
modes = modes(:,ordmode);

% Plot of mode shapes
scaleFactor = 2;
for ii = 1:3
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  
 

%% Damping Matrix

alfah = 1.3;
betah = 5e-7;

R = alfah*M + betah*K;

% Additional damper
rc = 50;
Rc = [rc -rc;-rc rc];
idfrc = [idf(7,2) idf(8,2)];
R(idfrc,idfrc)=R(idfrc,idfrc)+Rc;

RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in D
freq =0:0.05:75;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfDv = idf(7,2);
f0(idfDv)=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
idfCo = idf(5,1);
FRF_F_Co = xx(idfCo,:);

idfEv = idf(8,2);
idfDv = idf(7,2);
delta_l_DE = xx(idfDv,:)-xx(idfEv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_Co))
subplot(2,1,2)
plot(freq,angle(FRF_F_Co))

figure
subplot(2,1,1)
plot(freq,abs(delta_l_DE))
subplot(2,1,2)
plot(freq,angle(delta_l_DE))

%% Motion of Ax, effect on Cx and Ey

MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

xc = zeros(ndoc,1);
idfAo = idf(1,1)-ndof;
xc(idfAo) = 1;

clear A
for cc = 1:length(freq)
    A = -Om(cc)^2*MFF + sqrt(-1)*Om(cc)*RFF + KFF;
    B = -Om(cc)^2*MFC + sqrt(-1)*Om(cc)*RFC + KFC;
    xf(:,cc) = inv(A)*(-B*xc);
end

FRF_FA_Ev = xf(idfEv,:);
FRF_FA_Co = xf(idfCo,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_FA_Co))
subplot(2,1,2)
plot(freq,angle(FRF_FA_Co))

figure
subplot(2,1,1)
plot(freq,abs(FRF_FA_Ev))
subplot(2,1,2)
plot(freq,angle(FRF_FA_Ev))

%% Time response of Ev to Ao
T = 5;
Am = 0.1;
nHarm = 3;
t = 0.001:0.001:T;
N = length(t);

Amax = ones(length(t)/4,1).*3*Am/4;
Amin = ones(3*length(t)/4,1).*(-Am/4);
DisplA = [Amax; Amin];

figure
plot(t, DisplA)

% Fourier transform
fD = fft(DisplA,N);
df = 1/T;
nmax = N/2;
f = df*[0:(nmax-1)];
mod_fD(1) = abs(fD(1))/N;
mod_fD(2:nmax) = 2*abs(fD(2:nmax))/N;
phase_fD(1) = 0;
phase_fD(2:nmax) = angle(fD(2:nmax));

figure
subplot(211)
plot(f,mod_fD)
subplot(212)
plot(f,phase_fD)

% Frequency response
index_Harm = [2 3 4];
Om_Harm = 2*pi*f(index_Harm);
D_Harm = mod_fD(index_Harm);

d0 = zeros(ndoc,1);

for iHarm = 1:length(Om_Harm) % nHarm
    d0(idfAo)=D_Harm(iHarm);
    Afr=-Om_Harm(iHarm)^2*MFF + sqrt(-1)*Om_Harm(iHarm)*RFF + KFF;
    Bfr=-Om_Harm(iHarm)^2*MFC + sqrt(-1)*Om_Harm(iHarm)*RFC + KFC;
    x(:,iHarm) = inv(A)*-B*d0;
end

modx = abs(x);
phasex = angle(x);

% Time response
td = 0:0.001:10;
xtot = zeros(ndof,length(td));

for jHarm = 1:nHarm
    for it = 1:length(td)
        av = modx(:,jHarm);
        var = cos(Om_Harm(jHarm)*td(it)+phasex(:,jHarm));
        xtot(:,it) =  xtot(:,it) - av.*var;
    end
end

Ex = xtot(idfEv,:);
figure
plot(td,Ex)





