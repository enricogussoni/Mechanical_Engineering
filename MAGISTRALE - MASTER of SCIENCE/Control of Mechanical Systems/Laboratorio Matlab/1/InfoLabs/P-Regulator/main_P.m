clc
clear all
close all

%-----------------------------------------
% System data
 
g = 9.81;       %[m/s^2] gravity acceleration
 
m = 80;         %[kg] bar mass
J = 20;         %[kg*m^2] bar moment of inertia
L = 1;          %[m] bar length
 
k = 250;        %[N/m] stiffness coeff.
c = 25;         %[Ns/m] damping coeff.
 
% Generalized mass, damping and stiffness
ms = m*L^2 + J;
cs = 4*c*L^2;
ks = 4*k*L^2 - m*g*L;

num_Gm = [1];               %numerator
den_Gm = [ms cs ks];        %denominator
 
Gm = tf(num_Gm,den_Gm);     %TF of the passive system

disp('Passive system (no control)')
disp(' ')
Gm
disp('Poles of the passive system (no control)')
pole(Gm)


disp('Natural frequency of the passive system [Hz]:')
w0 = sqrt(ks/ms);
w0/2/pi
disp('Damping factor of the passive system [-]:')
h = cs/(2*ms*w0)


%-----------------------------------------
%-----------------------------------------
% Feedback control system

disp(' ')
disp('------------------------------------')
disp('Feedback control system')
disp(' ')
kp = input('Proportional gain: ')

% Open-loop TF
num_GH = [kp];              %numerator
den_GH = [ms cs ks];        %denominator
 
GH = tf(num_GH,den_GH)      %open-loop TF

disp('Poles of G(s)H(s):')
pole(GH)

%-------------------------
% Stability analysis - Undirect methods

% Nyquist criterion
figure(1)
hold on
nyquist(GH)

% Bode criterion
figure(2)
hold on
bode(GH)

% Phase and gain margin
figure
margin(GH)

%-------------------------
% Stability analysis - Direct methods

% Closed-loop TF
num_L = [kp];                   %numerator
den_L  = [ms  cs  (ks + kp)];	%denominator
 
L = tf(num_L,den_L)             %closed-loop TF

disp('Poles of L(s):')
pole(L)				

% Root locus
Gs = tf(num_GH/kp,den_GH);

figure
rlocus(Gs)

% Eigenvalues of the state matrix of the feedback control system
Ac = [-cs/ms  -(ks + kp)/ms
    1  0];
[Vc,Dc] = eig(Ac);
disp('Eigenvalues of the state matrix of the feedback control system [Ac]:')
diag(Dc)

%-------------------------
% Performance analysis

% Time domain analysis - Step response
figure(3)
hold on
step(L)

disp('Steady-state value of L(s):')
theta_inf = dcgain(L)

disp('Steady-state error:')
e_inf = 1 - theta_inf

% Frequency domain analysis - Bode diagrams of L(s)
figure(4)
hold on
bode(L)

disp('Bandwidth of the feedback control system [rad/s]:')
bandwidth(L)
