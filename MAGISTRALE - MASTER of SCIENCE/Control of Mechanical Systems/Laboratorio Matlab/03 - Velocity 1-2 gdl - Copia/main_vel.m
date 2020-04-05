clc
clear all
close all

%-----------------------------------------
% System data
 
Jm = 0.01;          % [kg*m2] motor inertia
Jr = 0.15;          % [kg*m^2] load inertia
r = 0.01;           % [Nms/rad] damping
kt = 1000;          % [Nm/rad] torsional spring

Cr = 20;            % [Nm] constant
% Cm is assumed to be Cm = (alfa*omega^2 + beta*omega + gamma)*y
Cm = [-0.01 1 10];  % [Nm] coefficients alfa,beta,gamma assuming Cm as a parabola for y=y0
y0 = 1;             % input variable at the equilibrium


% Equilibrium
omega = [0:0.01:110];
Cm_plot = Cm(1)*omega.^2 + Cm(2)*omega + Cm(3);
% figure
% plot(omega,Cm_plot)
% hold on
% plot(omega,Cr+2*r*omega,'r')

% Equilibrium values of omega
omega0 = roots(Cm-[0 2*r Cr]);
omega0 = sort(omega0);

disp('Steady-state values of omega [rad/s]:')
disp(num2str(omega0))

% Linearization
stability = 2;
switch stability
    case 1 %unstable case
        omega0 = omega0(1);
    case 2 %stable case
        omega0 = omega0(2);
end
a = 2*Cm(1)*omega0+Cm(2); % linearization with respect to omega
b = Cm(1)*omega0^2+Cm(2)*omega0+Cm(3); % linearization with respect to y

% equation of the system
n_dof = 2; % number of degrees of freedom
switch n_dof
    case 1
        M = [Jm+Jr];
        R = 2*r-a;
        K = 0;
        lambda = b;
    case 2
        M = [Jm 0; 0 Jr];
        R = [r-a 0; 0 r];
        K = [kt -kt; -kt kt];
        lambda = [b;0];
end

% Computation of the transfer functions: Method 1
s = tf('s');
system_pos1 = inv(M*s^2+R*s+K)*lambda;

% Computation of the transfer functions: Method 2 (using System state space)
A = [-inv(M)*R -inv(M)*K;eye(n_dof) zeros(n_dof)];
B = [inv(M)*lambda;zeros(n_dof,size(lambda,2))];
C = [zeros(n_dof) eye(n_dof)];
D = zeros(n_dof,size(lambda,2));
system_pos2 = tf(ss(A,B,C,D));

% Velocity output instead of position
system = minreal(system_pos2*s)

% Poles of the system
disp('Eigenvalues of A')
eig(A)
disp('Poles of system')
pole(system)


%-----------------------------------------
%-----------------------------------------
% Feedback control system

disp(' ')
disp('------------------------------------')
disp('Feedback control system')

% Open-loop TF
control_type = 2;
kp_vect = [0.1 0.5 1];

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
    
    feedback_variable = 1; % 1: co-located   2: non co-located
    
    Control = tf(num_C,den_C)
    G = series(Control,system);
    H = 1;

    GH = series(G(feedback_variable),H); %open-loop TF
    L = feedback(G,H,1,feedback_variable); %closed loop function
    
    disp('Open-loop function:')
    GH
    disp('Poles and zeros of G(s)H(s):')
    poles_GH = pole(GH)
    poles_GH_abs = abs(poles_GH)
    zero(GH)
    
    % Stability analysis from L
    disp('Closed-loop function:')
    L
    disp('Poles L(s):')
    pole(L)
    
    %-------------------------
    % Stability analysis

%     % Nyquist criterion
%     figure(100)
%     hold on
%     nyquist(GH)
% 
%     % Bode criterion: Phase and gain margin
%     figure(200)
%     hold on
%     margin(GH)
    
    %----------------------------------
    % Performance analysis
    % Response to a step reference
    figure(300)
    hold on
    step(L(feedback_variable),3)
    
    % Comparison between G (of the variable used as feedback) and L (H=1)
%     figure
%     bode(G(feedback_variable))
%     hold on
%     bode(L(feedback_variable))
%     title(['G and L for kp=' num2str(kp)])
%     
%     % Computation of steady-state error and bandwidth
%     theta_inf(count) = dcgain(L(feedback_variable));
%     e_inf(count) = 1 - theta_inf(count);
% 
%     % Frequency domain analysis - Bode diagrams of L(s)
%     bandw(count) = bandwidth(L(feedback_variable));
end

% figure(400)
% rlocus(GH/kp)

disp('Steady-state value of L(s) with respect to kp:')
disp(' ')
disp({'kp','Error','Bandw'})
disp(num2str([kp_vect' e_inf' bandw']))

