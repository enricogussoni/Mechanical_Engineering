%% Shink fit and shaft-hub connection

close all
clear all
clc

E = 200000;    %Nella FEM usa questo valore
v = 0.3;       % vm = va
UTS = 875;     % [MPa]
YS=545;

Cn = 2.15e6;   % [Nmm]
Y1 = 53487;    % [N]
f = 0.2;       % friction coefficient
Ks = 3;        % is a coefficient taking into account impacts and irregular working
               % conditions of the transmission
L = 177;       % [mm]
d = 212;       % [mm]
               
Clim = Cn*Ks;
Ylim = Y1*Ks;
               
%% From loads to required contact pressure
           
pc = 2*Clim/(f*pi*L*d^2);
tauc = pc*f;

py = Ylim/(f*pi*L*d);
tauy = py*f;

p = sqrt(tauc^2 + tauy^2)/f; % = pmin necessary to transmit the torque
                             % and to resist at the axial load
disp('Pressione minima:');
disp(p);

%% From geometric interference to available contact pressure

% dMmin = 212;     % [mm] H7
% dMmax = 212.046; % [mm] H7
% dAmin = 212.29;  % [mm] 212+0.29
% dAmax = 212.32;  % [mm] 212+0.32
% 
% imin = dAmin-dMmax;
% imax = dAmax-dMmin;

% Type 1

Dim = 212;
Dem = 420;
Dea = 212;
Dia = 60.5;

am = Dem/Dim;
aa = Dea/Dia;

% At the inner diameter of the hub
sigmar_m = -p;  
sigmat_m = p*(am^2+1)/(am^2-1);
epst_m1 = 1/E * (sigmat_m - v*sigmar_m);

% At the external diameter of the shaft
sigmar_a = -p;
sigmat_a = -p*(aa^2+1)/(aa^2-1);
epst_a1 = 1/E * (sigmat_a - v*sigmar_a);

imin_1 = (p*d/E)*((am^2+1)/(am^2-1)+(aa^2+1)/(aa^2-1));

% Type 2 (Grammel)

ri = 106;
r1 = 140;
r2 = 372.5;
re = 430;

h1 = 185.5;
h2 = 32;
h3 = 135;

xn=[0.7823 -20026 0.0358 -27905 0.1178 -21781];
x=xn*p;

sigmar_ri = x(1)+x(2)/ri^2;
sigmat_ri = x(1)-x(2)/ri^2;

epst_m2 = 1/E * (sigmat_ri - v*sigmar_ri);
epst_a2=epst_a1;

deltadm=epst_m2*d;
deltada=epst_a2*d;

imin_2=deltadm-deltada;
imin_2=p*d/E*(xn(1)*(1-v)-xn(2)/ri^2*(1+v)+(aa^2+1)/(aa^2-1)-v);

%% Effect of roughness on the interference
Rpm = 6.4e-3; %Rp=2*Ra
Rpa = 3.2e-3;

i_min1 = imin_1 + 2*Rpa + 2*Rpm;
i_min2 = imin_2 + 2*Rpa + 2*Rpm;
disp('Interferenza minima:');
disp(i_min1);
disp(i_min2);

%% comparison between standard and design
imin_dr=0.244;
imax_dr=0.32;

i_max_el = imax_dr -(2*Rpa + 2*Rpm);
i_min_el = imin_dr -(2*Rpa + 2*Rpm);

if i_min_el>=0.0009*d && i_max_el<=0.0015*d
    disp('Standard requirement is verified');
else
    disp('Minimum presssure does NOT verify standard');
end

%% Trovo effettiva pressione dato il design con massima inerferenza
tau_lim = YS;

%Type1
p1=i_max_el*E/d/((am^2+1)/(am^2-1)+(aa^2+1)/(aa^2-1));

sigmar_m = -p1;  
sigmat_m = p1*(am^2+1)/(am^2-1);

tau_max1 = max([abs(sigmar_m-sigmat_m),abs(sigmat_m),abs(sigmar_m)]);
eta1=tau_lim/tau_max1

%Type2
p2 = i_max_el*E/d/(xn(1)*(1-v)-xn(2)/ri^2*(1+v)+(aa^2+1)/(aa^2-1)-v);

x=xn*p2;
sigmar_ri = x(1)+x(2)/ri^2;
sigmat_ri = x(1)-x(2)/ri^2;

tau_max2 = max([abs(sigmar_ri-sigmat_ri),abs(sigmat_ri),abs(sigmar_ri)]);
eta2=tau_lim/tau_max2

disp('Interface contact pressure with maximum design interference:');
disp(p1);
disp(p2);

%% Determine the temperature variation that is needed to perform wheel installation
%  on the shaft.

deltadm = 0.64;    % [mm] = 2*imax
alpha = 0.000012;  % [1/°C]

epsilontheta = deltadm/d;

deltaT = epsilontheta/alpha

% Dalle slide, si può scaldare al massimo di 250 gradi, e a noi viene 251.
% Perciò facciamo che riscaldiamo di 250 e raffreddiamo l'albero.

%% Comparison with FEM results
Abaqus=[ 0.               104.125      
                2.95              76.2894     
                5.9               81.988      
                8.85              79.2686     
               11.8               56.7909     
               14.75              69.5908     
               17.7               72.6602     
               20.65              66.5886     
               23.6               68.3733     
               26.55              75.1285     
               29.5               72.8173     
               32.45              72.7331     
               35.4               77.3187     
               38.35              79.4453     
               41.3               78.9359     
               44.25              81.8186     
               47.2               85.0043     
               50.15              86.5053     
               53.1               88.3902     
               56.05              90.9068     
               59.                92.8809     
               61.95              94.9125     
               64.9               95.9181     
               67.85              97.9417     
               70.8              100.32       
               73.75             101.082      
               76.7              101.552      
               79.65             104.016      
               82.6              105.749      
               85.55             104.729      
               88.5              105.682      
               91.45             108.056      
               94.4              106.89       
               97.35             106.183      
              100.3              107.48       
              103.25             107.448      
              106.2              106.656      
              109.15             106.033      
              112.1              106.227      
              115.05             105.491      
              118.               104.365      
              120.95             103.216      
              123.9              101.516      
              126.85             100.96       
              129.8               98.9086     
              132.75              97.3508     
              135.7               94.9364     
              138.65              93.0385     
              141.6               91.8408     
              144.55              88.4617     
              147.5               86.2427     
              150.45              86.0047     
              153.4               83.5735     
              156.35              78.4848     
              159.3               79.9216     
              162.25              82.4169     
              165.2               72.142      
              168.15              75.3987     
              171.1               90.698      
              174.05             115.844      
              177               159.164          ];
ydistance=Abaqus(:,1);
pressure=Abaqus(:,2);
p1=@(x) p1+0*x;
p2=@(x) p2+0*x;
xdis=linspace(ydistance(1),ydistance(end),1000);

figure()
hold on
plot(ydistance,pressure);
plot(xdis,p1(xdis),'g');
plot(xdis,p2(xdis),'r');
legend('FEM results','Type1','Type2');
xlabel('True distance along path [mm]');
ylabel('Stress [MPa]')
title('Contact pressure');



