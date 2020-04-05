%% Esercizio MTU 2 DoF

%% Data
Jm = 0.01; %[Kg*M2]
J = 0.15;
r = 0.01; %[Ns/m], equivalente ressistance coeff.
kt = 1000; %[Nm/rad]
Tr = 20; %[Nm]
Tm = [-0.01 1 10]; % coefficienti della curva di coppia del motore = [alpha beta gamma]
% Tm = (alpha*omega^2 + beta*omega + gamma)*y

%% Equilibium
omega = [0:0.1:120];
Tm_plot = Tm(1)*omega.^2 + Tm(2)*omega + Tm(3);
Tr_plot = Tr+2*r*omega; % equivalent resistance torque

% figure(1)
% plot(omega,Tm_plot)
% hold on
% grid on
% plot(omega,Tr_plot,'r')
% hold off

omega0= roots([Tm(1) -2*r+Tm(2) Tm(3)-Tr]);
omega0=sort(omega0); % sort ordina in ordine crescente gli elementi del vettore

switch 2
    case 1
        omega0 = omega0(1); % unstable case
    case 2
        omega0 = omega0(2); % stable case
end

% Linearization of Tm
a = 2*Tm(1)*omega0+Tm(2); % cefficiente dello sviluppo in serie di Taylor
b = Tm(1)*omega0^2+Tm(2)*omega0+Tm(3);

%% Equation of motion of the system

switch 2
    case 1 % 1 dof
       n_dof = 1;
       M = Jm+J;
       R = 2*r-a;
       K = 0;
       lambda = b;
    case 2 % 2 dof
       n_dof = 2;
       M = [Jm 0 ; 0 J];
       R = [r-a 0 ; 0 r];
       K = [kt -kt ; -kt kt];
       lambda = [b ; 0];
end

% State space
A = [-inv(M)*R -inv(M)*K ; eye(n_dof) zeros(n_dof)];
autovaloriA = eig(A);
B = [inv(M)*lambda ; zeros(n_dof,1)];

% Laplace method 1
s = tf('s'); % s da ora in poi è identificata dal programma come la variabile di Laplace
system1_pos = inv(M*s^2+R*s+K)*lambda; % pos -> position

% Laplace method 2
% help ss
C = [zeros(n_dof) , eye(n_dof)]; % Misura della posizione (per coerenza del sistema)
D = zeros(n_dof,1); % zeros(number of measurements ; number of inputs)
system2_pos = tf(ss(A,B,C,D));

system_vel = minreal(series(system1_pos,s)); % minreal semplifica s nell'espressione del sistema (elimina zero e poli nello stesso punto)

%% Velocity control
kp = 1;
Ti = 1;
Td = 4;
feedback_variable = 2; % 1 => co-located
                       % 2 => non co-located

switch 2
    case 1 % Proportional
        control = kp;
    case 2 % Proportional Integrative
        control = kp*tf([Ti 1],[Ti 0]);
    case 3 % Proportional Derivative
         control = kp*tf([1/Td 1],[1/Td 0]); % Derivative is not so useful in this case
                                             % The derivative add one zero and
                                             % so the slope of the bode diagram
                                             % at high frequencies will be 0 =>
                                             % no attenuation
end

G = series(control,system_vel);
H = 1; % Ideal sensor
GH = series(G(feedback_variable,1),H); % espressione generale per 1 e 2 dof
L = feedback(G,H,1,feedback_variable); % espressione generale per 1 e 2 dof

% figure(2)
% margin(GH)

figure(3)
bode(G)
hold on
bode(L)

figure(4)
nyquist(GH)

figure(5)
rlocus(GH/kp)

figure(6)
step(L,1)
