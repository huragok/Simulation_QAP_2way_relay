clear all;
close all;
clc;

addpath('./functions');

%% 0. Load the data file that contains the test result
load('Test_2015730212237054') % 64QAM

%% 1. Simulation settings
N_batch = 5; % Number of batches,
N_per_batch = 1e4; % Number of monte-carlo run per batch, restricted by memory size
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

%% 3. Not let us test the ergodic mutual information
% EMI_analytical = cell(n_sigma2, 1);
EMI_MC = cell(n_sigma2, 1);

for i_sigma2 = 1 : n_sigma2
    tic

    % Compute the ergodic mutual information using our analytical approximated lower bound
%     EMI_analytical{i_sigma2} = zeros(M, 3);
%     EMI_analytical{i_sigma2}(:, 1) = get_EMI_lower_bound(constellation, map_noncore, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
%     EMI_analytical{i_sigma2}(:, 2) = get_EMI_lower_bound(constellation, map_core, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
%     EMI_analytical{i_sigma2}(:, 3) = get_EMI_lower_bound(constellation, map_QAP{i_sigma2}, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    
    
    % Compute the ergodic mutual information using Monte-Carlo simulation
    EMI_MC{i_sigma2} = zeros(M, 3);
    EMI_MC{i_sigma2}(:, 1) = get_EMI(constellation, map_noncore, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    EMI_MC{i_sigma2}(:, 2) = get_EMI(constellation, map_core, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    EMI_MC{i_sigma2}(:, 3) = get_EMI(constellation, map_QAP{i_sigma2}, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    
    toc;
    disp(['EMI simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
%     disp([' - EMI lower bounds, non-CoRe: ', num2str(EMI_analytical{i_sigma2}(:, 1)')])
%     disp([' - EMI lower bounds, CoRe: ', num2str(EMI_analytical{i_sigma2}(:, 2)')])
%     disp([' - EMI lower bounds, QAP: ', num2str(EMI_analytical{i_sigma2}(:, 3)')])
    disp([' - EMI emperical, non-CoRe: ', num2str(EMI_MC{i_sigma2}(:, 1)')])
    disp([' - EMI emperical, CoRe: ', num2str(EMI_MC{i_sigma2}(:, 2)')])
    disp([' - EMI emperical, QAP: ', num2str(EMI_MC{i_sigma2}(:, 3)')])
end

% EMI_analytical = reshape(cell2mat(EMI_analytical), M, 3 * n_sigma2);
EMI_MC = reshape(cell2mat(EMI_MC), M, 3 * n_sigma2);

%% Visualization
% The approximated EMI lowerbound
cmap = [0, 0, 0 ;0.5, 0, 1; 0, 0, 1; 1, 0, 0];
legend_item = cell(3 * M - 2, 1);
% h = figure;
% plot(dB_inv_sigma2, EMI_analytical(1, 1 : n_sigma2), 'k+-', 'linewidth', 2), hold on;
% legend_item{1} = 'TR0';
% for m = 2 : M
%     plot(dB_inv_sigma2, EMI_analytical(m, 1 : n_sigma2), '+-', 'Color', cmap(m, :), 'linewidth', 2);
%     plot(dB_inv_sigma2, EMI_analytical(m, n_sigma2 + 1 : 2 * n_sigma2), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
%     plot(dB_inv_sigma2, EMI_analytical(m, 2 * n_sigma2 + 1 : 3 * n_sigma2), 's-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;
%     
%     legend_item{3 * m - 4} = ['NM', num2str(m-1)];
%     legend_item{3 * m - 3} = ['CR', num2str(m-1)];
%     legend_item{3 * m - 2} = ['QAP', num2str(m-1)];
% end
% grid on;
% set(gca, 'Fontsize', 18);
% xlabel('1/\sigma^2(dB)'), ylabel('EMI');
% legend(legend_item, 'Location', 'northeast');
% saveas(h, ['EMI_noise_power_upperbound_new_', num2str(Q), 'QAM.fig']);

%% The empirical EMI
h = figure;
plot(dB_inv_sigma2, EMI_MC(1, 1 : n_sigma2), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    plot(dB_inv_sigma2, EMI_MC(m, 1 : n_sigma2), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    plot(dB_inv_sigma2, EMI_MC(m, n_sigma2 + 1 : 2 * n_sigma2), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    plot(dB_inv_sigma2, EMI_MC(m, 2 * n_sigma2 + 1 : 3 * n_sigma2), 's-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    
    legend_item{3 * m - 4} = ['NM', num2str(m-1)];
    legend_item{3 * m - 3} = ['CR', num2str(m-1)];
    legend_item{3 * m - 2} = ['QAP', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('EMI');
legend(legend_item, 'Location', 'northeast');
saveas(h, ['EMI_noise_power_MonteCarlo_new_', num2str(Q), 'QAM.fig']);