clc
clear all
close all

%-----------------------------------------
% System data
 
g = 9.81;       %[m/s^2] gravity acceleration
 
m = 80;         %[kg] bar mass
J = 20;         %[kg*m^2] bar moment of inertia
l = 1;          %[m] bar length

stability = 1;
switch stability
    case 1 % stable
        k = 250;        %[N/m] stiffness coeff.
    case 2 % unstable (statically)
        k = 125;
    case 3 % ks=0
        k = 196.2;
end

c = 25;         %[Ns/m] damping coeff.
 
% Generalized mass, damping and stiffness
ms = m*l^2 + J;
cs = 4*c*l^2;
ks = 4*k*l^2 - m*g*l;

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

% Open-loop TF
control_type = 2;
kp_vect = [100 300 1000 5000];

count = 0;
for kp=kp_vect
    count = count+1;
    switch control_type
        case 0 % P
            
            num_C = kp;
            den_C = 1;
        case 1 % PD
            Td = 4;

            num_C = [Td*kp kp];
            den_C = [1];
        case 2 % PI
            Ti = 0.9;
            
            num_C = [Ti*kp kp];
            den_C = [Ti 0];
        otherwise

    end
    
    Control = tf(num_C,den_C)
    G = series(Control,Gm);
    H = 1;

    GH = series(G,H); %open-loop TF
    L = feedback(G,H); %closed loop function
    
    disp('Open-loop function:')
    GH
    disp('Poles and zeros of G(s)H(s):')
    poles_GH = pole(GH)
    poles_GH_abs = abs(poles_GH)
    zero(GH)
    
    % Stability analysis from L
    disp('Closed-loop function:')
    L
    disp('Poles and zeros of L(s):')
    pole(L)
    zero(L)
    
    %-------------------------
    % Stability analysis

    % Nyquist criterion
    figure(100)
    hold on
    nyquist(GH)

    % Bode criterion: Phase and gain margin
    figure(200)
    hold on
    margin(GH)
    
    %----------------------------------
    % Performance analysis
    % Response to a step reference
    figure(300)
    hold on
    step(L)
    
    % Comparison between G and L (H=1)
    figure
    bode(G)
    hold on
    bode(L)
    title(['G and L for kp=' num2str(kp)])
    
    % Computation of steady-state error and bandwidth
    theta_inf(count) = dcgain(L);
    e_inf(count) = 1 - theta_inf(count);

    % Frequency domain analysis - Bode diagrams of L(s)
    bandw(count) = bandwidth(L);
end

figure(400)
rlocus(GH/kp)

disp('Steady-state value of L(s) with respect to kp:')
disp(' ')
disp({'kp','Error','Bandw'})
disp(num2str([kp_vect' e_inf' bandw']))

