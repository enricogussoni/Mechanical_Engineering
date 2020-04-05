clear
close all
clc

%% Structure properties
mp = 4;
EAp = 4e7;
EJp = 4e3;
mr = 1;
EAr = 1e5;
EJr = 1e3;

fmax = 20;
SF = 2;
ommax = 2*pi*fmax*SF;

Lmaxp = sqrt(pi^2/ommax * sqrt(EJp/mp))
Lmaxr = sqrt(pi^2/ommax * sqrt(EJr/mr))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('040915');


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

%% Partitioning of "FF" Mass and Stiffness Matrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes, Om2] = eig(MFF\KFF);
om2 = diag(Om2);
frq = sqrt(om2)/2/pi;
om = sqrt(om2);

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

h1 = 0.02;
h2 = 0.03;
h = [h1 h2];

% coeff = [1/(2*om(1)) om(1)/2 ; 1/(2*om(2)) om(2)/2];
frq_damp=frqord(1:2);
coeffs = [1./(2*2*pi*frq_damp) (2*pi*frq_damp/2)];
ab = coeffs\h';

alfah = ab(1);
betah = ab(2);

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in k
freq = 0:0.01:20 ;
Om = 2*pi*freq;

f0 = zeros(ndof,1);
idfk = idf(6,2);
f0(idfk)=-1;

for ii=1:length(freq)    
%     A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
%     xx(:,ii)=A\f0;
xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF)*f0;
end
idfyj = idf(7,2);
FRF_Fk_yj = xx(idfyj,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_Fk_yj))
subplot(2,1,2)
plot(freq,angle(FRF_Fk_yj))

%% Reaction force in C due to which results from a periodic vertical force 
% applied in node k consisting in the superposition of two harmonic components

MCF = M(ndof+1:end,1:ndof);
KCF = K(ndof+1:end,1:ndof);
RCF = R(ndof+1:end,1:ndof);

om1 = 3;
om2 = 4;
A1 = 500;
A2 = 500;
phi1 = 0;
phi2 = pi/4;

f01 = zeros(ndof,1);
f02 = f01;

f01(idfk)=500;
f02(idfk)=500;

for cc=1:length(freq)
    xf1(:,cc) = (-Om(cc)^2*MFF +sqrt(-1)*Om(cc)*RFF + KFF)^-1 * f01;
    xf2(:,cc) = (-Om(cc)^2*MFF +sqrt(-1)*Om(cc)*RFF + KFF)^-1 * f02;
    R1(:,cc) = (-Om(cc)^2*MCF +sqrt(-1)*Om(cc)*RCF + KCF)*xf1(:,cc);
    R2(:,cc) = (-Om(cc)^2*MCF +sqrt(-1)*Om(cc)*RCF + KCF)*xf2(:,cc);
    R(:,cc) = R1(:,cc)*exp(sqrt(-1)*(om1*cc+phi1)) + R2(:,cc)*exp(sqrt(-1)*(om2*cc+phi2));
end

% for ii=1:length(freq)    
%     A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
%     B = -Om(ii)^2*MCF+sqrt(-1)*Om(ii)*CCF+KCF;
%     F0(idof_y6) =500*cos(3*ii)+500*cos(4*ii+45) ;
%     rr(:,ii)=(B)*inv(A)*F0;
% end

idfCv = idf(1,2)-ndof;
FRF_Fk12_Cv = R(idfCv,:);

figure
subplot(211)
plot(freq,abs(FRF_Fk12_Cv))
subplot(212)
plot(freq,angle(FRF_Fk12_Cv))

% Time response
T = 10;
Fmax = 500; %[N]
% n_harmonics = ;
fsample=1000;
dt = 1/fsample;
t = dt:dt:T;
N = length(t);

Force = 500*cos(om1*t + phi1)+500*cos(om2*t + phi2);
figure
plot(t,Force);grid
xlabel('Time [s]');ylabel('Periodic Force [N]')

% Fourier analysis with fft
fForce = fft(Force);
nmax = N/2;
delta_f = 1/T;
freq = 0:delta_f:delta_f*(nmax-1);
mod_Force(1) = fForce(1)/N;
mod_Force(2:nmax) = 2*abs(fForce(2:nmax))/N;
phase_Force(1) = 0;
phase_Force(2:nmax) = angle(fForce(2:nmax));

figure
subplot(211)
bar(freq,mod_Force)
subplot(212)
plot(freq,phase_Force,'o')

% Frequency response 
index_harm = [1:10]; 
freq_harm = freq(index_harm);
om_harm = 2*pi*freq_harm;
Force_harm = mod_Force(index_harm);

F0 = zeros(ndof,1);

% superposition of the effects of the n harmonics
clear xf rr
for iHarm = 1:length(index_harm)   %nharmonics
    F0(idfCv) = Force_harm(iHarm);
    xf(:,iHarm) = inv(-om_harm(iHarm)^2*MFF +1i*om_harm(iHarm)*RFF +KFF)^-1 *F0; 
    rr(:,iHarm) = (-om_harm(iHarm)^2*MCF +1i*om_harm(iHarm)*RCF +KCF)*xf(:,iHarm);
end

modr = abs(rr);
phaser = angle(rr);

% Time domain
td = [0:0.01:T];
r = zeros(ndoc,length(td));

for iHarm = 1:length(index_harm)
    for it = 1:length(td)
        r(:,it) = r(:,it)+modr(:,iHarm).*cos(om_harm(iHarm)*td(it)+phaser(:,iHarm));
        % A * cos(om*t + phi)
    end
end

rC = r(idfCv,:);
figure
plot(td,rC)

%% Displacement of A due to static and distributed vertical load

p = -200;
% gammaAB = gamma(10);
gammaAB = atan(1/4);
pn = p*cos(gammaAB);
pt = p*sin(gammaAB);
p0 = zeros(ndof,1);

idfr = [idf(8,3) idf(11,3)];
idfv = idf(8:11,2);
idfo = idf(8:11,1);

vertical_r = pn.*l(7:10)/2 * cos(gammaAB) + pt.*l(7:10) * sin(gammaAB);
orizontal_r = pn.*l(7:10)/2 * sin(gammaAB) + pt.*l(7:10) * cos(gammaAB);
flectional_r = pn^2.*[l(7) l(10)]/12;

p0(idfr)=flectional_r;
p0(idfv)=vertical_r;
p0(idfo)=orizontal_r;

xx = KFF\p0;

idfAo = idf(8,1);
dispAo = xx(idfAo,1)
idfAv = idf(8,2);
dispAv = xx(idfAv,1)

