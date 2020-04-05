%% Data

L=1; % [m]
m=80; % [kg]
J=2; % [kg*m2]
r=25; %[Ns/m]
g=9.81; %[m/s2]
k=196.2; %[N/m]

% from indicatore1=1:2


% Equivalent phisical properties

mx=J+m*L^2;
rx=4*r*L^2;
kx=4*k*L^2 - m*g*L;
lambda= L;

%% System

roots([mx rx kx]); % trova le radici del polinomio risolvente
A=[-rx/mx -kx/mx; 1 0];
eig(A);

system = tf(lambda,[mx rx kx]); % transfer function -> Continuous-time transfer function are defined in Laplace space

P=pole(system);
[Z,gain]=zero(system);

% Control block
ti_vet=1:2:10;

for ii=1:length(td_vet)
    kp=500;
    td=0.1;
    ti=ti_vet(ii);

    switch 3 % metere 1,2 o 3 a seconda del caso da studiare
        case 1 % Proportional
            control = kp;
        case 2 % Proportional-Derivative
            control = kp*tf([td 1],1);
        case 3 % Proportional-Integral
            control = kp*tf([ti 1],[ti 0]); %the constant term is 0
    end

    G=series(control,system); % equivalent TF
    H=1; % ideal sensor
    GH=series(G,H); % open-loop TF
    closeL= feedback(G,H); % closed-loop TF

    %% Diagrams

    figure(1)
    %bode(GH)
    margin(GH) % for phase margin and gain margin
    hold on
    % bode(closeL)

    figure(2)
    nyquist(GH)
    hold on % per vedere il plot al variare di kp in kp_vett

    figure(3)
    rlocus(GH/kp) % root locus varing the gain Kp that has to be escluded from the function
                  % per maggiori output leggi l'help
    hold on

    figure(4) % Step-function response
    step(closeL,20)
    hold on
    
end