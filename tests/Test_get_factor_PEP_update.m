% Script to test the functionality of get_factor_PEP_update function and
% compare its result with the true factor in the form of E[exp()] using
% Monte-Carlo simulation

clear all;
close all;
clc;

addpath('../functions');

%% 1. Simulation settings
% Constellation specification
Nbps = 4;
type_mod = 'QAM';
pwr = 1; 

% Node S, R, D power, channel power and noise power specification
% For now we assume that:
%   S and D are using the same unit power 
%   all the 3 node S, R, D's AWGN noise have the same power sigma2
% We assume that the S and D are separated by distance 1. The variance of
% the relay channel is d ^ -nu where nu is the pathloss factor
dB_inv_sigma2 = 0; % 1/sigma2 in dB
Pr = 2; % Power at the relay
d1 = 0.5; % Distance between S and R
d2 = 0.5; % Distance between R and D
nu = 3; % Pathloss factor

N = 50000; % Number of monte-carlo run

%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, pwr);

sigma_sqr = 10 ^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr / (beta_sr + beta_rd + sigma_sqr_r)); % The power normalization factor

gamma_sr = abs(sqrt(beta_sr / 2) * (randn(N, 1) + 1i * randn(N, 1))) .^ 2; % Generate the channel absolute value square
gamma_rd = abs(sqrt(beta_rd / 2) * (randn(N, 1) + 1i * randn(N, 1))) .^ 2;
%% 3. Test some  E
dist_sqr = [0, 0.2, 0.5, 1, 2, 4, 8, 16];
n_dist_sqr = length(dist_sqr);

E_analytical = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r);

E_montecarlo = zeros(1, n_dist_sqr);
for i_dist_sqr = 1 : n_dist_sqr
    E_montecarlo(i_dist_sqr) = mean(exp(-g ^ 2 / 4 * dist_sqr(i_dist_sqr) * gamma_sr .* gamma_rd ./ (sigma_sqr_d + g ^ 2 * sigma_sqr_r * gamma_rd)));
end

E_accurate_montecarlo = zeros(1, n_dist_sqr);
for i_dist_sqr = 1 : n_dist_sqr
    E_accurate_montecarlo(i_dist_sqr) = mean(qfunc(sqrt(g ^ 2 / 2 * dist_sqr(i_dist_sqr) * gamma_sr .* gamma_rd ./ (sigma_sqr_d + g ^ 2 * sigma_sqr_r * gamma_rd))));
end

%% Visualization
h = figure;
plot(dist_sqr, E_analytical / 2, 'bs-', 'linewidth', 2), hold on;
plot(dist_sqr, E_montecarlo / 2, 'ro--', 'linewidth', 2);
plot(dist_sqr, E_accurate_montecarlo, 'k-', 'linewidth', 2);
grid on;
set(gca, 'fontsize', 16);
xlabel('\epsilon'), ylabel('E');
legend('Chernoff analytical', 'Chernoff Monte-Carlo', 'Monte-Carlo');
saveas(h, 'Test_get_factor_PEP_update.fig')
