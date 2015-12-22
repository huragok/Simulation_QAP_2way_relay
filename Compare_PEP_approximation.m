clear all;
close all;
clc;

addpath('./functions');


%% 1. Simulation settings
% Constellation specification
Nbps = 4; % 4, 5, 6
type_mod = 'QAM';

% Node S, R, D power, channel power and noise power specification
% For now we assume that:
%   S and D are using the same unit power 
%   all the 3 node S, R, D's AWGN noise have the same power sigma2
% We assume that the S and D are separated by distance 1. The variance of
% the relay channel is d ^ -nu where nu is the pathloss factor
% We also assume that the channel and noise are stationary across
% transmissions

%dB_inv_sigma2 = [1 : 1.5 : 13]; % 1/sigma2 in dB, 64 QAM
dB_inv_sigma2 = [-2 : 1.5 : 10]; % 1/sigma2 in dB, 64 QAM
p_Pr = 0.5; % this portion of the total power of 4 is allocated to the relay. The rest are divided eqaully between the 2 end nodes
d = [0.5, 0.5]; % Distance between S and R, R and D

nu = 3; % Pathloss factor
% We set M = 2, since for lar Number of retransmission

N_batch = 5; % Number of batches,
N_per_batch = 1e6; % Number of monte-carlo run per batch, restricted by memory size
seed = 8;

% 64QAM
% p0 = 0; q0 = 6;
% p1 = 0; q1 = 48;

% 16QAM
p0 = 0; q0 = 3;
p1 = 0; q1 = 15;

%% 2. Initialization: generate and save all test cases
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes are the same n_sigma2-by-1 vector
n_sigma2 = length(dB_inv_sigma2);

Pr = 4 * p_Pr; % Power at the relay
Pt = 2 * (1 - p_Pr); % Power at each of the 2 transmitter

beta_sr = d(:, 1) .^ -nu; % n_d-by-1 vector
beta_rd = d(:, 2) .^ -nu;

dist_sqr0 = Pt * abs(constellation(p0 + 1) - constellation(q0 + 1)) ^ 2;
dist_sqr1 = Pt * abs(constellation(p1 + 1) - constellation(q1 + 1)) ^ 2;

pep_num = zeros(n_sigma2, 1); % PEP evaluated according to its definition in (7) numerically
pep_chernoff_approx = zeros(n_sigma2, 1);  % PEP evaluated using Chernoff bound and the approximation by Proposition 1
pep_doubleexp_approx = zeros(n_sigma2, 1); % PEP evaluated using more accurate sum of two exponential approximation and the approximation by Proposition 1
pep_chernoff_num = zeros(n_sigma2, 1); % PEP evaluated using Chernoff bound with E evaluated numerically
pep_doubleexp_num = zeros(n_sigma2, 1); % PEP evaluated using the more accurate sum of two exponential approximation with E evaluated numerically
rng(seed);

for i_batch = 1 : N_batch
    % generate the random channels for the first transmission
    delta1_0 = abs(sqrt(beta_sr / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    gamma2_0 = abs(sqrt(beta_rd / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    
    % generate the random channels for the first retransmission
    delta1_1 = abs(sqrt(beta_sr / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    gamma2_1 = abs(sqrt(beta_rd / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
   
    for i_sigma2 = 1 : n_sigma2
        g = sqrt(Pr / (beta_sr * Pt + beta_rd * Pt + sigma_sqr(i_sigma2))); % The power normalization factor
        
        alpha_0 = sqrt(Pr ./ (delta1_0 * Pt + gamma2_0 * Pt + sigma_sqr(i_sigma2))); % The power normalization factor
        sigma_sqr_tilde_0 = sigma_sqr(i_sigma2) + sigma_sqr(i_sigma2) * alpha_0 .^ 2 .* gamma2_0; % The effective noise power at the destination
        tmp0 = alpha_0 .^ 2 .* gamma2_0 .* delta1_0./ (sigma_sqr_tilde_0 .^ 2);

        alpha_1 = sqrt(Pr ./ (delta1_1 * Pt + gamma2_1 * Pt + sigma_sqr(i_sigma2))); % The power normalization factor
        sigma_sqr_tilde_1 = sigma_sqr(i_sigma2) + sigma_sqr(i_sigma2) * alpha_1 .^ 2 .* gamma2_1; % The effective noise power at the destination
        tmp1 = alpha_1 .^ 2 .* gamma2_1 .* delta1_1./ (sigma_sqr_tilde_1 .^ 2);
        
        pep_num(i_sigma2) = pep_num(i_sigma2) + mean(qfunc(sqrt(dist_sqr0 * tmp0 / 2 + dist_sqr1 * tmp1 / 2)));
        pep_chernoff_num(i_sigma2) = pep_chernoff_num(i_sigma2) + 0.5 * mean(exp(-dist_sqr0 * tmp0 / 4)) * mean(exp(-dist_sqr1 * tmp1 / 4));
        pep_doubleexp_num(i_sigma2) = pep_doubleexp_num(i_sigma2) + mean(exp(-dist_sqr0 * tmp0 / 4)) * mean(exp(-dist_sqr1 * tmp1 / 4)) / 12 + mean(exp(-dist_sqr0 * tmp0 / 3)) * mean(exp(-dist_sqr1 * tmp1 / 3)) / 4;
        
        pep_chernoff_approx(i_sigma2) =...
        pep_chernoff_approx(i_sigma2) +...
        0.5 * get_factor_PEP_update_param(dist_sqr0, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 4)...
            * get_factor_PEP_update_param(dist_sqr1, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 4);
        
        pep_doubleexp_approx(i_sigma2) =...
        pep_doubleexp_approx(i_sigma2) +...
        1 / 12 * get_factor_PEP_update_param(dist_sqr0, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 4)...
               * get_factor_PEP_update_param(dist_sqr1, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 4)...
        +1 / 4 * get_factor_PEP_update_param(dist_sqr0, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 3)...
               * get_factor_PEP_update_param(dist_sqr1, beta_sr, beta_rd, g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), 3);
    end
end

pep_num = pep_num / N_batch;
pep_chernoff_approx = pep_chernoff_approx / N_batch;
pep_doubleexp_approx = pep_doubleexp_approx / N_batch;
pep_chernoff_num = pep_chernoff_num / N_batch;
pep_doubleexp_num = pep_doubleexp_num / N_batch;

figure;
semilogy(dB_inv_sigma2, pep_num, 'k+:', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2, pep_doubleexp_num, 'rs-.', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_chernoff_num, 'ro-.', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_doubleexp_approx, 'bs-', 'linewidth', 2);
semilogy(dB_inv_sigma2, pep_chernoff_approx, 'bo-', 'linewidth', 2);
legend('Qfunc+Num', 'ExpSum+Num', 'Chernoff+Num', 'ExpSum+Approx', 'Chernoff+Approx');
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('PEP');