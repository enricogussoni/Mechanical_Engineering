%% Generic material and geometry data
clc
clear

% Evaluation of Kt

Kt = 1.16;
sigma_u=600;  % = Rm
sigma_e=200;  % Numerator of sigma_e_first (from standard)
rho = 0.3^2;  % from grafic
r=15;         % [mm], radiud of the notch (from design)
q = 1/(1+sqrt(rho/r));
Kf = 1+q*(Kt-1);
sigma_e_prime = 200/Kf;  % 200 Nmm from standard

%% Calculation of the nominal stress with De Saint Venant beam theory

%Mx   = 60245e3;   % [Nmm], NON E' QUESTO; VA VALUTATO CON SOLO PESO
Mtot = 60385e3;   % [Nmm], " " "
Mt   =  -2142e3;   % [Nmm], " " "
Y1   = 53487;     % [N],   " " "

D = 180;          % [mm], external (in the smaller section)
d = 60.5;         % [mm], internal (worst case)
Jp = pi/32 * (D^4 - d^4);

sigmay_sf = 26.71; % più o meno

sigmay_Mx   = 93.53;
sigmay_Mtot = 32*Mtot*D/pi/(D^4 - d^4);
sigmay_ax   = Y1/(pi/4 * (D^2 - d^2)); % Non c'è invece da valutare Kft e
                                         % Kfn perché il torcente e
                                         % l'azione assiale sono costanti.
sigmay_braking = sigmay_ax + sigmay_sf;
tauxy       = 16*Mt*D/pi/(D^4 - d^4);

%% Global method (it considers Kt explicitly)
% Analytical (Or from graphs)

Sigma_a_m = [0      0      0;
             0   sigmay_Mx 0;
             0      0      0];  % sigma versor for alternate stress in motion

Sigma_m_m = [0      0      0;
             0   sigmay_braking 0;
             0      0      0];  % sigma versor for mean stress in motion
         
Sigma_a_b = [0       0       0;
             0   sigmay_Mtot 0;
             0       0       0]; % sigma versor for alternate stress while braking
         
Sigma_m_b = [0       tauxy         0;
             tauxy sigmay_braking 0;
             0       0             0]; % sigma versor for mean stress while braking

% The effect of mean stress is considered trough Haig diagrams and Goodman criterion
% In motion
Im = sigmay_braking;
sigma_s=sigmay_Mx;
eta_mot=sigma_e_prime*(1-Im/sigma_u)/sigma_s

% While braking
Im = sigmay_braking;
sigma_s = sigmay_Mtot;
eta_br = sigma_e_prime*(1-Im/sigma_u)/sigma_s
