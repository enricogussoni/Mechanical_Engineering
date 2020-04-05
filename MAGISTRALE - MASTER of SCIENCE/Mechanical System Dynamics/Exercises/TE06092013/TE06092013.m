clear
close all
clc

%% Structure properties
EA = 7.4e7;
EJ = 2e4;
m = 2.8;

fmax = 25;
SC = 2;
ommax = 2*pi*SC*fmax;

Lmax = sqrt(pi^2/ommax * sqrt(EJ/m))


%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE06092013');


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

% Additional springs
k = 2000;
kadd = [k -k;-k k];

% idf1o = idf(1,1);
% idf2o = idf(2,1);
% idf3o = idf(3,1);
% K([idf1o idf2o],[idf1o idf2o])=K([idf1o idf2o],[idf1o idf2o])+kadd;
% K([idf2o idf3o],[idf2o idf3o])=K([idf2o idf3o],[idf2o idf3o])+kadd;

idf12o = idf([1 2],1);
idf23o = idf([2 3],1);
K(idf12o,idf12o)=K(idf12o,idf12o)+kadd;
K(idf12o,idf12o)=K(idf12o,idf12o)+kadd;

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

alfah = 0.6;
betah = 4e-5;

R = alfah*M + betah*K;
RFF = R(1:ndof,1:ndof);


%% Frequency Response Function

% Force applied in A
freq = 0:0.01:25;
Om = 2*pi*freq;

f0 = zeros(ndof,1) ;
idfAo = idf(11,1);
f0(idfAo)=1;

clear ii
for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f0;
end
FRF_F_A = xx(idfAo,:);
idfBo = idf(13,1);
FRF_F_B = xx(idfBo,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F_A))
subplot(2,1,2)
plot(freq,abs(FRF_F_B))

%% Calcolare la risposta in frequenza da 0-25 Hz con passo 0.01 Hz dello 
% spostamento verticale del punto A per uno spostamento orizzontale dei 
% vincoli C e D di ampiezza unitaria e fra loro in contro-fase.

MFC=M(1:ndof,ndof+1:end);
RFC=R(1:ndof,ndof+1:end);
KFC=K(1:ndof,ndof+1:end);

freq = 0:0.01:25;
Om = 2*pi*freq;

xc = zeros(ndoc,1);
idfDo = idf(7,1)-ndof;
idfCo = idf(5,1)-ndof;
xc(idfDo)=1;
xc(idfCo)=-1;

for cc=1:length(freq)
    A = -Om(cc)^2*MFF + sqrt(-1)*Om(cc)*RFF + KFF;
    B = -Om(cc)^2*MFC + sqrt(-1)*Om(cc)*RFC + KFC;
    xf(:,cc) = inv(A)*-B*xc;
end

FRF_CDo_Ao = xf(idfAo,:);

figure
subplot(211)
plot(freq,abs(FRF_CDo_Ao))
subplot(212)
plot(freq,angle(FRF_CDo_Ao))

%% Calcolare il valore (costante) della reazione verticale nei carrelli 
% O1, O2, O3 prodotto dal peso proprio della struttura nella posizione di equilibrio.
m = 2.8;
p = m*9.81;
p0 = zeros(ndof,1);
l_or = l(9);
l_vert = l(1)*2;

% idf_vert_pos = idf([5 7 9 10 12 15],2);
% p0(idf_vert_pos)=ones(length(idf_vert_pos))*p*l_or/2;
% p0(idf(14,2))=p*l_or;
% idf_vert_neg = idf([3 5 7 9 19 12 14 15 1 2],2);
% p0(
% p0(idf(1,2))=-p*l_vert/2;
% p0(idf(2,2))=-p*l_vert;
% p0(idf(3,1))=-p*l_vert/2;
p0(idf(15,2))=p*l_or/2 - p*l_vert/2;
p0(idf(15,3))=p*l_or^2/12;
p0(idf(14,2))=p*l_or - p*l_vert/2;
p0(idf(5,2))=p*l_or/2;
p0(idf(5,3))=p*l_or^2/12;
p0(idf(12,2))=p*l_or/2 - p*l_vert/2;
p0(idf(12,3))=p*l_or^2/12;
p0(idf(7,2))=p*l_or/2 - p*l_vert/2;
p0(idf(7,3))=p*l_or^2/12;
p0(idf(9,2))=p*l_or/2 - p*l_vert/2;
p0(idf(9,3))=p*l_or^2/12;
p0(idf(10,2))=p*l_or/2 - p*l_vert/2;
p0(idf(9,3))=p*l_or^2/12;

MCF=M(ndof+1:end,1:ndof);
RCF=M(ndof+1:end,1:ndof);
KCF=M(ndof+1:end,1:ndof);

A = KFF;
xf(:,1)=A\p0;
B = KCF;
rr(:,1)=B*xf(:,1)

rrO1 = rr(idf(1,1))


