%% Exam Simulation
% 6. Assuming that the external torque C(t) has the periodic form reported in Fig.4,
% with a period T = 3s and Cmax = 500 Nm, compute the steady-state response of
% the horizontal displacement A, neglecting any contribution of the external load for
% harmonics higher than the third.

%% Periodic external torque

T = 3;
Cmax = 500;
fsample = 1000;
dt = 1/fsample;
t = dt:dt:T;
N = length(t);
Torque = zeros(1,N);
Torque(1:end/2) = Cmax;
Torque(end/2+1:end) = -Cmax;

figure
plot(t,Torque);
title('One period of the signal')
% All the odd armonic will be equal to 0

fTorque = fft(Torque);
nmax = N/2;
delta_f = 1/T;
freq = delta_f*[0:(nmax-1)];
Torque_mod(1) = fTorque(1)/N;  % Torque modulus at 0 Hz
Torque_mod(2:end) = 2*abs(fTorque(2:end))/N;
Torque_phase(1)=0;
Torque_phase(2:end)= angle(fTorque(2:end));



