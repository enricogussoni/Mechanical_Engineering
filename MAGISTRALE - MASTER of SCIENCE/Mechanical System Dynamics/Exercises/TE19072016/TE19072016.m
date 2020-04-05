clear
close all
clc

%% Structure properties 
E = 2.06e11;  % [N/m^2] -> Pa (non MPA)
rho = 7800;% [kg/m^3]
Ab = 3.912e-3 ;  % [m^2]
Jb = 3.892e-5;
Ar = 0.01155;
Jr = 4.82e-4 ;  % [m^4]
EJb = E*Jb; % [Nm^2]
EJr = E*Jr;
EAb = E*Ab; % [N]
EAr = E*Ar; 
mr=rho*Ar; % [kg/m]
mb=rho*Ab;

fmax = 20; 
SC = 1.5;
Ommax = 2*pi*SC*fmax;

Lmaxr = sqrt(pi^2/Ommax * sqrt(EJr/mr))
Lmaxb = sqrt(pi^2/Ommax * sqrt(EJb/mb))

%% Load Structure Data

[file_name,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam] = MeccFEM2_loadstructure('TE19072016');


%% Plot undeformed structure

figure
MeccFEM2_plotStructure(position,l,gamma,xy)
xlabel('x [m]'); ylabel('y [m]')


%% Check IDB and ndof

MeccFEM2_DoFsTable(idf)


%% Assembly of Mass and Stiffness Matrices

ndof_total = 3*nnod;
[M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_total);
ndoc=size(M,1)-ndof;


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
scaleFactor = 6;
for ii = 1:5
    mode = modes(:,ii);
    figure
    MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(ii) ': Freq [Hz]=' num2str(frqord(ii))])
end  

%% Damping Matrix
alfah = 0.8;
betah = 3.0e-5;

R = alfah*M + betah*K;

RFF = R(1:ndof,1:ndof);

%% Frequency Response Function

% Force applied in A (F1) -> response in Cv
freq = 0:0.01:20;
Om = 2*pi*freq;

f1 = zeros(ndof,1);
idfAv = idf(19,2);
f1(idfAv)=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f1;
end

idfCv = idf(9,2);
FRF_F1_Cv = xx(idfCv,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F1_Cv))
subplot(2,1,2)
plot(freq,angle(FRF_F1_Cv))

% Force applied in B (F2) -> response in Co
f2 = zeros(ndof,1);
idfBo = idf(4,1);
f1(idfBo)=1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*RFF+KFF;
    xx(:,ii)=A\f1;
end

idfCo = idf(9,1);
FRF_F1_Co = xx(idfCo,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_F1_Co))
subplot(2,1,2)
plot(freq,angle(FRF_F1_Co))

%% Response to a displacement of 1 and 7 (oriz. and in phase)
% -> displacement in Co

% Partitioning of "FC" Matrices
MFC = M(1:ndof,ndof+1:end);
RFC = R(1:ndof,ndof+1:end);
KFC = K(1:ndof,ndof+1:end);

xc = zeros(ndoc,1); 

idf1o=idf(1,1)-ndof; 
xc(idf1o) = 1;
idf17o=idf(17,1)-ndof;
xc(idf17o) = 1;

for ii=1:length(freq)    
    A = -Om(ii)^2*MFF  +1i*Om(ii)*RFF  +KFF;
    B = -Om(ii)^2*MFC  +1i*Om(ii)*RFC  +KFC;
    xf(:,ii)= inv(A)*-B*xc;
end

idfCo = idf(9,1);
FRF_eq_Co = xf(idfCo,:);

figure
subplot(2,1,1)
plot(freq,abs(FRF_eq_Co))
subplot(2,1,2)
plot(freq,angle(FRF_eq_Co))

%% Compute the static vertical displacement of C due to a snow load equal 
% to 9400 N/m applied on the top beams.

p0 = 9400;
Floc = zeros(ndof,1);
Fglo = zeros(ndof,1);

gl = gamma(6);
lambdal = [cos(gl) sin(gl) 0; -sin(gl) cos(gl) 0; 0 0 1];
pol = p0*sin(gl);
pvl = p0*cos(gl);

Floc(idf(5,:),1) = [pol*l(5)/2 pvl*l(5)/2 pvl*l(5)^2/12];
Fglo(idf(5,:),1) = Floc(idf(5,:),1)'*lambdal;

Floc(idf(6,:),1) = [pol*l(6)/2 pvl*l(6)/2 0];
Fglo(idf(6,:),1) = Floc(idf(6,:),1)'*lambdal;

Floc(idf(7,:),1) = [pol*l(7)/2 pvl*l(7)/2 0];
Fglo(idf(7,:),1) = Floc(idf(7,:),1)'*lambdal;

Floc(idf(8,:),1) = [pol*l(8)/2 pvl*l(8)/2 0];
Fglo(idf(8,:),1) = Floc(idf(8,:),1)'*lambdal;

Flocl(idf(9,:),1) = [pol*l(8)/2 pvl*l(8)/2 pvl*l(8)^2/12];
Fglo(idf(9,:),1) = Flocl(idf(9,:),1)'*lambdal;

gr = gamma(10);
lambdar = [cos(gr) sin(gr) 0; -sin(gr) cos(gr) 0; 0 0 1];
por = -p0*sin(gr);
pvr = p0*cos(gr);

Flocr(idf(9,:),1) = [por*l(9)/2 -pvr*l(9)/2 -pvr*l(9)/12];
Fglo(idf(9,:),1) = Fglo(idf(9,:),1) + (Flocr(idf(9,:),1)'*lambdal)'; 

Flocr(idf(10,:),1) = [por*l(10)/2 -pvr*l(10)/2 -pvr*l(10)/12];
Fglo(idf(10,:),1) = Floc(idf(10,:),1)'*lambdar;

Flocr(idf(11,:),1) = [por*l(11)/2 -pvr*l(11)/2 -pvr*l(11)/12];
Fglo(idf(11,:),1) = Floc(idf(11,:),1)'*lambdar;

Flocr(idf(12,:),1) = [por*l(12)/2 -pvr*l(12)/2 -pvr*l(12)/12];
Fglo(idf(12,:),1) = Floc(idf(12,:),1)'*lambdar;

Flocr(idf(13,:),1) = [por*l(13)/2 -pvr*l(13)/2 -pvr*l(13)/12];
Fglo(idf(13,:),1) = Floc(idf(13,:),1)'*lambdar;

    
    A = KFF;
    xxn = A\Fglo;

Displ_snow_Cv = xxn(idfCv);

% p = 9400;
% % g1 = gamma(5)
% % g2 = gamma(15)
% 
% g = gamma(6);
% po = -p*sin(g);
% pv = p*cos(g);
% 
% Floc = zeros(ndof,1);
% Fglo = zeros(ndof,1);
% 
% Floc(idb(6,:),1) = [po*l(5)/2 -pv*l(5)/2 -pv*l(5)^2/12];
% lambda = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];
% Fglo(idb(6,:),1) = Floc(idb(6,:),1)'*lambda;
% 
% Floc(idb(7,:),1) = [po*(l(5)+l(6))/2, -pv*(l(5)+l(6))/2, pv*(l(5)^2-l(6)^2)/12];
% Fglo(idb(7,:),1) = Floc(idb(7,:),1)'*lambda;
% 
% Floc(idb(8,:),1) = [po*(l(6)+l(7))/2, -pv*(l(6)+l(7))/2, pv*(l(6)^2-l(7)^2)/12];
% Fglo(idb(8,:),1) = Floc(idb(8,:),1)'*lambda;
% 
% Floc(idb(9,:),1) = [po*(l(7)+l(8))/2, -pv*(l(7)+l(8))/2, pv*(l(7)^2-l(8)^2)/12];
% Fglo(idb(9,:),1) = Floc(idb(9,:),1)'*lambda;
% 
% Floc(idb(10,:),1) = [po*l(8)/2, -pv*l(8)/2, pv*l(8)^2/12];
% Fglo(idb(10,:),1) = Floc(idb(10,:),1)'*lambda;
% 
% g = gamma(17);
% 
% Floc(idb(6,:),1) = [po*l(14)/2, pv*l(14)/2, pv*l(14)^2/12];
% lambda = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];
% Fglo(idb(6,:),1) = Fglo(idb(6,:),1)' + Floc(idb(6,:),1)'*lambda;
% 
% Floc(idb(15,:),1) = [po*(l(14)+l(15))/2, pv*(l(14)+l(15))/2, pv*(l(15)^2-l(14)^2)/12];
% Fglo(idb(15,:),1) = Floc(idb(15,:),1)'*lambda;
% 
% Floc(idb(16,:),1) = [po*(l(15)+l(16))/2, pv*(l(15)+l(16))/2, pv*(l(16)^2-l(15)^2)/12];
% Fglo(idb(16,:),1) = Floc(idb(16,:),1)'*lambda;
% 
% Floc(idb(17,:),1) = [po*(l(16)+l(17))/2, pv*(l(16)+l(17))/2, pv*(l(17)^2-l(16)^2)/12];
% Fglo(idb(17,:),1) = Floc(idb(17,:),1)'*lambda;
% 
% Floc(idb(18,:),1) = [po*l(17)/2, pv*l(17)/2, pv*l(17)^2/12];
% Fglo(idb(18,:),1) = Floc(idb(18,:),1)'*lambda;
% 
% xstat = inv(KFF)*(Fglo);
% vertdispC = xstat(idb(6,2))
% hordispC = xstat(idb(6,1))


