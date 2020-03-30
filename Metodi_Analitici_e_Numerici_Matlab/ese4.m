clear all
close all
clc

%% PUNTO 1

% Dati sperimentali
sigma   = [0.18   0.3   0.5    0.6    0.72  0.75   0.8   0.9   1.0   ];
epsilon = [0.0005 0.001 0.0013 0.0015 0.002 0.0045 0.006 0.007 0.0085];

figure(1)
% axes('FontSize',12)
plot(sigma, epsilon, 'ko', 'linewidth', 2)
title('Dati sperimentali')
xlabel('Sforzo')
ylabel('Deformazione')

%% Interpolazione di Lagrange (IL)
sigma_dis  = linspace(min(sigma), max(sigma), 1000);
grado      = length(sigma) - 1;
PL         = polyfit(sigma, epsilon, grado); 
epsilon_IL = polyval(PL, sigma_dis);

figure(2)
axes('FontSize',12)
plot(sigma, epsilon, 'ko', sigma_dis, epsilon_IL, 'r', 'linewidth', 2)
xlabel('Sforzo')
ylabel('Deformazione')
title('Dati sperimentali & Interp. Pol. Lagrange')
   
%% Interpolazione composita lineare  (ICL)
epsilon_ICL = interp1(sigma, epsilon, sigma_dis);

figure(3)
axes('FontSize',12)
plot(sigma, epsilon, 'ko', sigma_dis, epsilon_ICL, 'r', 'linewidth', 2)
xlabel('Sforzo')
ylabel('Deformazione')
title('Dati sperimentali & Interp. Comp. Lineare')
   

%% PUNTO 2

sigma_v        = [0.4 0.65];
epsilon_IL_v   = polyval(PL, sigma_v);
epsilon_ICL_v  = interp1(sigma, epsilon, sigma_v);

% Print
fprintf('Sigma:         %f   %f\n', sigma_v)
fprintf('Epsilon IL:   %f  %f\n', epsilon_IL_v)
fprintf('Epsilon ICL:   %f   %f\n', epsilon_ICL_v)
