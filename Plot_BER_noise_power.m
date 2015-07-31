clear all;
close all;
clc;

addpath('./functions');

%% 0. Load the data file that contains the test result
%load('./data/Test_201578104535979.MAT') % 16QAM
%load('Test_201571310354702') % 32QAM
load('Test_201573020554265') % 64QAM


%% 1. Simulation settings
N_batch = 5; % Number of batches,
N_per_batch = 1e6; % Number of monte-carlo run per batch, restricted by memory size
seed = 8;

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;

p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
d1 = test_cases(1).param_origin.d1; % Distance between S and R
d2 = test_cases(1).param_origin.d2; % Distance between R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

n_sigma2 = length(test_cases);
dB_inv_sigma2 = zeros(1, n_sigma2);
for i_sigma2 = 1 : n_sigma2
    dB_inv_sigma2(i_sigma2) = test_cases(i_sigma2).param_origin.dB_inv_sigma2;
end

%% 2. Initialization
Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor

map_noncore = get_map_noncore(Q, M);
map_core = get_map_core(Q, M);
map_QAP = cell(n_sigma2, 1);
for i_sigma2 = 1 : n_sigma2
    map_QAP{i_sigma2} = [1 : Q; test_cases(i_sigma2).map];
end

%% 3. Not let us test the bit error rate
BER_analytical = cell(n_sigma2, 1);
BER_MC = cell(n_sigma2, 1);

for i_sigma2 = 1 : n_sigma2
    tic

    % Compute the bit error rate using our analytical upper bound
    BER_analytical{i_sigma2} = zeros(M, 3);
    BER_analytical{i_sigma2}(:, 1) = get_BER_upper_bound(constellation, map_noncore, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    BER_analytical{i_sigma2}(:, 2) = get_BER_upper_bound(constellation, map_core, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    BER_analytical{i_sigma2}(:, 3) = get_BER_upper_bound(constellation, map_QAP{i_sigma2}, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    
    % Compute the bit error rate using Monte-Carlo simulation
    BER_MC{i_sigma2} = zeros(M, 3);
    BER_MC{i_sigma2}(:, 1) = get_BER(constellation, map_noncore, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC{i_sigma2}(:, 2) = get_BER(constellation, map_core, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC{i_sigma2}(:, 3) = get_BER(constellation, map_QAP{i_sigma2}, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    
    toc;
    disp(['BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - BER upper bounds, non-CoRe: ', num2str(BER_analytical{i_sigma2}(:, 1)')])
    disp([' - BER upper bounds, CoRe: ', num2str(BER_analytical{i_sigma2}(:, 2)')])
    disp([' - BER upper bounds, QAP: ', num2str(BER_analytical{i_sigma2}(:, 3)')])
    disp([' - BER emperical, non-CoRe: ', num2str(BER_MC{i_sigma2}(:, 1)')])
    disp([' - BER emperical, CoRe: ', num2str(BER_MC{i_sigma2}(:, 2)')])
    disp([' - BER emperical, QAP: ', num2str(BER_MC{i_sigma2}(:, 3)')])
end

BER_analytical = reshape(cell2mat(BER_analytical), M, 3 * n_sigma2);
BER_MC = reshape(cell2mat(BER_MC), M, 3 * n_sigma2);

save(['BER_noise_power_', num2str(Q), 'QAM.mat'], 'dB_inv_sigma2', 'BER_analytical', 'BER_MC');
%% Visualization
% The BER upperbound
cmap = [0, 0, 0 ;0.5, 0, 1; 0, 0, 1; 1, 0, 0];
legend_item = cell(3 * M - 2, 1);
h = figure;
semilogy(dB_inv_sigma2, BER_analytical(1, 1 : n_sigma2), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    semilogy(dB_inv_sigma2, BER_analytical(m, 1 : n_sigma2), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    semilogy(dB_inv_sigma2, BER_analytical(m, n_sigma2 + 1 : 2 * n_sigma2), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical(m, 2 * n_sigma2 + 1 : 3 * n_sigma2), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 4} = ['NM', num2str(m-1)];
    legend_item{3 * m - 3} = ['CR', num2str(m-1)];
    legend_item{3 * m - 2} = ['QAP', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
ylim([1e-6, 10]), xlim([8, 28])
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_noise_power_upperbound_', num2str(Q), 'QAM.fig']);

%% The empirical BER
h = figure;
semilogy(dB_inv_sigma2, BER_MC(1, 1 : n_sigma2), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    semilogy(dB_inv_sigma2, BER_MC(m, 1 : n_sigma2), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    semilogy(dB_inv_sigma2, BER_MC(m, n_sigma2 + 1 : 2 * n_sigma2), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_MC(m, 2 * n_sigma2 + 1 : 3 * n_sigma2), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 4} = ['NM', num2str(m-1)];
    legend_item{3 * m - 3} = ['CR', num2str(m-1)];
    legend_item{3 * m - 2} = ['QAP', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
ylim([1e-6, 10]), xlim([8, 28])
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_noise_power_MonteCarlo_', num2str(Q), 'QAM.fig']);