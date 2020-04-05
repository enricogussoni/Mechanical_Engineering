clear all
close all
clc



%% Structure properties


m1  = 4;    %[kg/m]
m2  = 1;    %[kg/m]

M1 = 4;        %[kg]
M2 = 50;        %[kg]

J1 = 0.04;      %[kgm^2]
J2 = 0.1;       %[kgm^2]

EA1 = 4e7;   %[N]
EA2 = 1e5;   %[N]
EJ1 = 4e3;   %[Nm^2]
EJ2 = 1e3;   %[Nm^2]

% Maximum length of finite element
fmax=20;
coef=2;
ommax=coef*(fmax*2*pi);
Lmax_fe_1 = sqrt(pi^2/ommax*sqrt(EJ1/m1))
Lmax_fe_2 = sqrt(pi^2/ommax*sqrt(EJ2/m2))


%% Load Structure Data
close all
[file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure('040915');

% Plot undeformed structure
figure
MeccFEM2_plotStructure(position,xy,m,EA,EJ);
xlabel('x [m]'); ylabel('y [m]')

% Check idf and ndof
MeccFEM2_DoFsTable(idf)
fprintf('Number of DoFs: %i\n',ndof)


%% Assembly of Mass and Stiffness Matrices

ndof_tot = 3*nnod;
[M,K]=MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof_tot);

%% Contribution due to k1
% idof_k1=idf(3,1);
% K(idof_k1,idof_k1) = K(idof_k1,idof_k1)+k1;
% 
%% Contribution due to k2
% K_k2 = [ k2 -k2 ;
%         -k2  k2 ];
%     
% idof_k2=idf([5,6],2);
% K(idof_k2,idof_k2)=K(idof_k2,idof_k2)+K_k2;
% 
%% Contribution due to k3
% K_k3_local = [ k3 0 0 -k3 0 0;
%                0  0 0  0  0 0;
%                0  0 0  0  0 0;
%               -k3 0 0  k3 0 0;
%                0  0 0  0  0 0;
%                0  0 0  0  0 0;];
% alpha_k3 = 3/4*pi;
% lambda_k3 = [ cos(alpha_k3) sin(alpha_k3) 0 ; 
%              -sin(alpha_k3) cos(alpha_k3) 0 ;
%               0      0      1 ];
%           
% Lambda_k3 = [lambda_k3  zeros(3,3) ; 
%              zeros(3,3) lambda_k3  ];
%          
% K_k3_global = Lambda_k3'*K_k3_local*Lambda_k3;
% 
% idof_n2 = idf(2,:);
% idof_n4 = idf(4,:);
% idof_k3 = [idof_n2 idof_n4];
% 
% K(idof_k3,idof_k3)=K(idof_k3,idof_k3)+K_k3_global;

%% Contribution due to lumped masses
% idof_n4 = idf(4,:);
% idof_n5 = idf(5,:);
% idof_n6 = idf(6,:);
% MM1 = diag([M1 M1 J1]);
% MM2 = diag([M2 M2 J2]);
% 
% M(idof_n4,idof_n4) = M(idof_n4,idof_n4)+MM1;
% M(idof_n5,idof_n5) = M(idof_n5,idof_n5)+MM2;
% M(idof_n6,idof_n6) = M(idof_n6,idof_n6)+MM1;

% Partitioning of "FF" Mass and Stiffness Matrices
MFF=M(1:ndof,1:ndof);
KFF=K(1:ndof,1:ndof);


%% Analysis of Natural Frequencies and Mode Shapes

[modes,Om2]=eig(MFF\KFF);
om2=diag(Om2);
frq=sqrt(om2)/(2*pi);

% Sort in ascending order frequencies and mode shapes
[frq,ordmode]=sort(frq);
modes = modes(:,ordmode);

nMod = sum(frq<fmax);

% Plot of mode shapes
fscale=1;
for i=1:3
    mode=modes(:,i);
    figure;
    MeccFEM2_plotDeformedStructure(mode,fscale,incidence,l,gamma,position,idf,xy);
    xlabel('x [m]'); ylabel('y [m]')
    title(['Mode ', num2str(i) ': Freq [Hz]=' num2str(frq(i))])
end


%% Modal mass and modal stiffness

Phi = modes(:,1:3);
Mmodal_3f = Phi'*MFF*Phi; % for the first 3 mode shapes
Kmodal_3f = Phi'*KFF*Phi; % for the first 3 mode shapes

% fprintf('Modal Mass      (1st to 3rd mode)  : [%d %d %d]\n', diag(Mmodal_3f));
% fprintf('Modal Stiffness (1st to 3rd mode)  : [%d %d %d]\n', diag(Kmodal_3f));
% 
% freq_calc = sqrt(diag(Kmodal_3f)./diag(Mmodal_3f))./(2*pi);
% fprintf('Natural frequency (1st to 3rd mode)  : [%i %i %i]\n', frq(1:3));
% fprintf('Natural frequency with modal approach (1st to 3rd mode)  : [%i %i %i]\n', freq_calc);

%% Damping Matrix
% Structural isteretic damping
h1 = 0.02;
h2 = 0.03;

h=[h1;h2];

frq_damp=frq(1:2);
coeffs = [1./(2*2*pi*frq_damp) (2*pi*frq_damp/2)];

alpha_beta = (coeffs'*coeffs)^-1*coeffs'*[h1;h2];
alpha_beta = coeffs\h;

alpha = alpha_beta(1);
beta = alpha_beta(2);

C=alpha*M+beta*K;

% Contribution due to c2
% C_c2 = [ c2 -c2 ; 
%         -c2 c2  ];
%     
% idof_c2 = idf([5,6],2);
% C(idof_c2,idof_c2)=C(idof_c2,idof_c2)+C_c2;

% Extract the C_FF matrix
CFF=C(1:ndof,1:ndof);

% Plot of h=h(f)
h_computed=alpha./(2*2*pi*frq_damp) + beta.*(2*pi*frq_damp)/2;

ffs=linspace(0,fmax,1000);
h_fun=alpha./(2*2*pi*ffs) + beta.*(2*pi*ffs)/2;

figure
plot(ffs,h_fun,'b',frq_damp,h_computed,'ro')
title('Non-dimensional damping ratio')
xlabel('Freq. [Hz]')
grid on


alpha = alpha_beta(1);
beta = alpha_beta(2);
C=alpha*M+beta*K;

CFF=C(1:ndof,1:ndof);


%% Frequency Response Function

freq=[0:0.1:fmax].';
Om=2*pi*freq;
F0 = zeros(ndof,1);
idof_k_vert = idf(6,2);

F0(idof_k_vert) = -1;

for ii=1:length(freq)    
    xx(:,ii)=inv(-Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF)*F0;
end

idof_y_vert = idf(7,2);
FRF_YV_F = xx(idof_y_vert,:);


figure
subplot(211)
plot(freq,abs(FRF_YV_F))
grid on
title('FRF: Vertical Displacement of C vs Force')
ylabel(['|Y_C/F| [m/N]'])
subplot(212)
plot(freq,angle(FRF_YV_F))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

%% Frequency Response Function - Distributed force
%Compute the static vertical and horizontal displacements of point A under the constant distributed 
%vertical load p=200 N/m that is acting only on the upper inclined beam and write the computed values
%in the table at the back of this paper.

f0 = zeros(ndof,1);
idof_y8 = idf(8,2);
idof_theta8 = idf(8,3);
idof_y7 = idf(7,2);
idof_theta7 = idf(7,3);
idof_y6 = idf(6,2);
idof_theta6 = idf(6,3);

idof_y9 = idf(9,2);
idof_theta9 = idf(9,3);
idof_y10 = idf(10,2);
idof_theta10 = idf(10,3);
idof_y11 = idf(11,2);
idof_theta11 = idf(11,3);

l_fe1 = (0.25^2+1)^0.5;
p0 = 200;

f0(idof_y8,1) = -p0*l_fe1/2;
f0(idof_theta8,1) = -p0*l_fe1^2/12;
f0(idof_y7,1) = -p0*l_fe1/2 *2;
f0(idof_theta7,1) = 0; % =-p0*l_fe^2/12 + p0*l_fe^2/12

l_fe2 = ((0.5^2+4)^0.5)/3;

f0(idof_y6,1) = -p0*l_fe1/2+-p0*l_fe2/2;
f0(idof_theta6,1) = 0;

l_fe2 = ((0.5^2+4)^0.5)/3;

f0(idof_y9,1) = -p0*l_fe2/2 *2;
f0(idof_theta9,1) = 0; % =-p0*l_fe^2/12 + p0*l_fe^2/12

f0(idof_y10,1) = -p0*l_fe2/2 *2;
f0(idof_theta7,1) = 0; % =-p0*l_fe^2/12 + p0*l_fe^2/12

f0(idof_y11,1) = -p0*l_fe2/2;
f0(idof_theta11,1) = p0*l_fe2^2/12;

for ii=1:length(freq)   
    A = -Om(ii)^2*MFF  +1i*Om(ii)*CFF   +KFF;
    xf(:,ii)=A\f0;
end


idof_XA = idf(8,1);
FRF_XA_p = xf(idof_XA,:);

idof_YA = idf(8,2);
FRF_YA_p = xf(idof_YA,:);

figure
subplot(211)
plot(freq,abs(FRF_XA_p))
grid on
title('FRF: Horizontal Displacement of A vs Distributed force')
ylabel(['|Y_A/p_0| [m/N/m]'])
subplot(212)
plot(freq,angle(FRF_XA_p))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on

figure
subplot(211)
plot(freq,abs(FRF_YA_p))
grid on
title('FRF: Vertical Displacement of A vs Distributed force')
ylabel(['|Y_A/p_0| [m/N/m]'])
subplot(212)
plot(freq,angle(FRF_YA_p))
ylabel(['\Psi [rad]'])
xlabel('Freq [Hz]')
grid on


%% Compute the constraint vertical force in point C which results from a periodic vertical 
%force applied in node k consisting in the superposition of two harmonic components: ?1=3Hz (A1=500N, ?1=0deg)
%and ?2=4Hz (A2=500N, ?2=45deg). Plot the time history (10 s) of the computed steady state vertical force 






% MCF = M(ndof+1:end,1:ndof);
% CCF = C(ndof+1:end,1:ndof);
% KCF = K(ndof+1:end,1:ndof);
% freq=[0:0.1:fmax].';
% Om=2*pi*freq;
% F0 = zeros(ndof,1);
% idof_y6 = idf(6,2);
% 
% 
% for ii=1:length(freq)    
%     A = -Om(ii)^2*MFF+sqrt(-1)*Om(ii)*CFF+KFF;
%     B = -Om(ii)^2*MCF+sqrt(-1)*Om(ii)*CCF+KCF;
%     F0(idof_y6) =500*cos(3*ii)+500*cos(4*ii+45) ;
%     rr(:,ii)=(B)*inv(A)*F0;
% end
% 
% idof_vert_c = idf(1,2) -ndof;
% 
% FRF_VF_F = rr(idof_vert_c,:);
% 
% figure
% subplot(211)
% plot(freq,abs(FRF_VF_F))
% grid on
% title('FRF: Vertical Constraint Force in C vs Force')
% ylabel(['|V_F/F| [m/N]'])
% subplot(212)
% plot(freq,angle(FRF_VF_F))
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')
% grid on
% 
% 
% % %% Periodic External force
% % 
% % T = 10;                %[s]
% % Fmax = 500;             %[N]
% % fsample1 = 3;         %[Hz]
% % dt1 = 1/fsample1;         % Time step
% % fsample2 = 4;         %[Hz]
% % dt2 = 1/fsample2;         % Time step
% % time_vec1 = [dt1:dt1:T];   % Time vector
% % time_vec2 = [dt2:dt2:T];
% % N1 = length(time_vec1);
% % N2 = length(time_vec2);
% % 
% % 
% % 
% % Force = ;
% % figure
% % plot(time_vec1,Force);grid
% % xlabel('Time [s]');ylabel('Periodic Force [N]')
% % axis([0 time_vec1(end) -Fmax*1.05 Fmax*1.05])
% % 
% % %Fourier Analysis with fft
% % Force_f = fft(Force,N);
% % nmax = N/2;
% % delta_freq = 1/T;                   %Frequency step 
% % freq_vec = delta_freq*[0:(nmax-1)]; %Frequency vector
% % Force_mod(1) = Force_f(1)/N;      %Torque modulus at 0 Hz 
% % Force_mod(2:nmax)=2/N*abs(Force_f(2:nmax));   %Torque modulus for frequencies > 0 Hz
% % Force_phase(1) = 0;
% % Force_phase(2:nmax) = angle(Force_f(2:nmax));          %Torque phase in radiant
% % 
% % figure
% % subplot(211)
% % bar(freq_vec,Force_mod,0.1);
% % set(gca,'Xlim',[0 500]);
% % subplot(212)
% % plot(freq_vec,Force_phase,'o')
% % 
% % %Frequency response
% % index_3_contributions = [2 4 6]; %even harmonic contributions and mean value are 0!
% % Ome_vec = 2*pi*freq_vec(index_3_contributions);
% % Force_f_3c = Force_mod(index_3_contributions).*exp(sqrt(-1)*Force_phase(index_3_contributions));
% % 
% % F0 = zeros(ndof,1);
% % idof_XA = idf(4,1);
% % for ii=1:length(Ome_vec)
% %     F0(idof_XA) = Force_f_3c(ii);
% %     xxp(:,ii)=inv(-Ome_vec(ii)^2*MFF+sqrt(-1)*Ome_vec(ii)*CFF+KFF)*F0;
% % end
% % 
% % % Time domain
% % t = [0:0.001:2];
% % xtot = zeros(ndof,length(t));
% % 
% % for jjj = 1:length(Ome_vec)
% %     for iii=1:length(t)
% %         xtot(:,iii)=xtot(:,iii)-Ome_vec(jjj)^2*abs(xxp(:,jjj)).*cos(Ome_vec(jjj)*t(iii)+angle(xxp(:,jjj)));
% %     end
% % end
% % idof_XB = idf(6,1);
% % figure
% % plot(t,xtot(idof_XB,:))
% % grid
% % xlabel('[s]')
% % ylabel('$\ddot{X}_B  [m/s^2]$','interpreter','latex')
% % title('Steady-state response of $\ddot{X}_B$','interpreter','latex')
% % 
% % %
% % dfu_dx=1/l(10); % 10 è l'elemento finito di interesse (FE5_1.inp!!!!)
% % idof_n3 = idf(3,:);
% % idof_n11 = idf(11,:);
% % lambda = [cos(gamma(10)) sin(gamma(10)) 0; -sin(gamma(10)) cos(gamma(10)) 0; 0 0 1];
% % X3L = lambda*xx(idof_n3,:);
% % X11L = lambda*xx(idof_n11,:);
% % load_ass_n3 = EA2*dfu_dx*(X3L(1,:)-X11L(1,:)); %K*Dl
% % 
% % %oppure (considero la K come EA/Ltot
% % dfu_dx=1/(2*l(10)); % 10 è l'elemento finito di interesse (FE5_1.inp!!!!)
% % 
% % load_ass_n3_opzione = EA2*dfu_dx*X3L(1,:); 
% % 
% % figure
% % subplot(211)
% % plot(freq,abs(load_ass_n3),freq,abs(load_ass_n3_opzione),'r')
% % grid on
% % title('FRF: Assial load in D vs Force')
% % ylabel(['|Ax Load/F| [N/N]'])
% % subplot(212)
% % plot(freq,angle(load_ass_n3),freq,angle(load_ass_n3_opzione),'r')
% % ylabel(['\Psi [rad]'])
% % xlabel('Freq [Hz]')
% % grid on
% 
% 
