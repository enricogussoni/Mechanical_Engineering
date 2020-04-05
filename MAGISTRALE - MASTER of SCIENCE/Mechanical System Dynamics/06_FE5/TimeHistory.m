%% Periodic External Torque
% Compute the steady-state response of the horizontal acceleration of B due 
% to a periodic force applied in the horizontal direction in A. The time 
% history of the force is reported in Figure 2 (T = 0.1 s and Fmax = 1E5 N). 
% Neglect any contribution of the external force for harmonics higher than 
% the third.


T = 0.1;
Fmax = 1e5;
nHarm = 3;
fsamp=1000;
dt = 1/fsamp;
t = dt:dt:T;
N = length(t);

Force = [interp1([0,T/4],[0,Fmax],dt:dt:T/4), ...
         interp1([T/4, 3/4*T],[Fmax,-Fmax],T/4+dt:dt:3/4*T), ...
         interp1([3*T/4, T],[-Fmax, 0],3/4*T+dt:dt:T)];

figure 
plot(t,Force)

%Fourier Analysis with fft
fForce = fft(Force,N);
nmax = N/2;
fres = 1/T;
f = 0:fres:fres*(nmax-1);
mod_fForce = abs(fForce);
phase_fForce = angle(fForce);
mod_fForce(1) = mod_fForce(1)./N;
mod_fForce(2:nmax) = 2.*mod_fForce(2:nmax)./N;
phase_fForce(1) = 0;

figure
subplot(211)
bar(f,mod_fForce(1:nmax))
subplot(212)
plot(f,phase_fForce(1:nmax),'o')

Force_mod = mod_fForce(1:nmax);
Force_phase = phase_fForce(1:nmax);

%Frequency response
index_Harm = [2:2:2*nHarm];
f_Harm = f(index_Harm);
Om_Harm = 2*pi*f_Harm;
Force_Harm = Force_mod(index_Harm).*exp(sqrt(-1)*Force_phase(index_Harm));

F0 = zeros(ndof,1);
idfAo = idf(4,1);

for ii = 1:length(f_Harm) % = length(Om-Harm) -> scorro le frequenze
    F0(idfAo)=Force_Harm(ii); % sovrapposizione degli effetti per la 
                              % forzante ad ogni frequenza considerata
    A = -Om_Harm(ii)^2*MFF + sqrt(-1)*Om_Harm(ii)*RFF + KFF;
    xx(:,ii) = A\F0;
end

% Time domain
td = 0:0.001:10;
xt = zeros(ndof,length(td)); % matrice delle accelerazioni: time history di ogni gdl

for it=1:length(td)
    for iHarm=1:length(Om_Harm)
        mean = -Om_Harm(iHarm)^2*abs(xx(:,iHarm));
        variation = cos(Om_Harm(iHarm)*td(it)+angle(xx(:,iHarm)));
        xt(:,it)=xt(:,it)-mean.*variation; % x = A*cos(om*t+phi)
    end
end

idf_Bo = idf(6,1);
figure
plot(t,xtot(idof_XB,:))