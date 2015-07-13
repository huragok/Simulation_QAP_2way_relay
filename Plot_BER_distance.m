clear all;
close all;
clc;

addpath('./functions');

%% 0. Load the data file that contains the test result
load('./data/Test_20157810591291.MAT') % 16QAM
%load('./data/Test_201551519134039.MAT') % 32QAM
%load('./data/Test_201551519847716') % 64QAM

%% 1. Simulation settings
N_batch = 5; % Number of batches,
N_per_batch = 1e6; % Number of monte-carlo run per batch, restricted by memory size
seed = 8;

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
dB_inv_sigma2 = test_cases(1).param_origin.dB_inv_sigma2;
p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

n_d = length(test_cases);
d1 = zeros(n_d, 1);
d2 = zeros(n_d, 1);
for i_d = 1 : n_d
    d1(i_d) = test_cases(i_d).param_origin.d1;
    d2(i_d) = test_cases(i_d).param_origin.d2;
end

%% 2. Initialization
Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 .^ -nu;
beta_rd = d2 .^ -nu;

g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor

map_noncore = get_map_noncore(Q, M);
map_seddik = get_map_seddik(Q, M);
map_QAP = cell(n_d, 1);
for i_d = 1 : n_d
    map_QAP{i_d} = [1 : Q; test_cases(i_d).map];
end

%% 3. Not let us test the bit error rate
BER_analytical = cell(n_d, 1);
BER_MC = cell(n_d, 1);

for i_d = 1 : n_d
    tic

    % Compute the bit error rate using our analytical upper bound
    BER_analytical{i_d} = zeros(M, 3);
    BER_analytical{i_d}(:, 1) = get_BER_upper_bound(constellation, map_noncore, beta_sr(i_d), beta_rd(i_d), g((i_d)), sigma_sqr_d, sigma_sqr_r);
    BER_analytical{i_d}(:, 2) = get_BER_upper_bound(constellation, map_seddik, beta_sr(i_d), beta_rd(i_d), g((i_d)), sigma_sqr_d, sigma_sqr_r);
    BER_analytical{i_d}(:, 3) = get_BER_upper_bound(constellation, map_QAP{i_d}, beta_sr(i_d), beta_rd(i_d), g((i_d)), sigma_sqr_d, sigma_sqr_r);
    
    % Compute the bit error rate using Monte-Carlo simulation
    BER_MC{i_d} = zeros(M, 3);
    BER_MC{i_d}(:, 1) = get_BER(constellation, map_noncore, beta_sr(i_d), beta_rd(i_d), Pr, Pt, Pt, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed);
    BER_MC{i_d}(:, 2) = get_BER(constellation, map_seddik, beta_sr(i_d), beta_rd(i_d), Pr, Pt, Pt, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed);
    BER_MC{i_d}(:, 3) = get_BER(constellation, map_QAP{i_d}, beta_sr(i_d), beta_rd(i_d), Pr, Pt, Pt, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed);
    
    toc;
    disp(['Distance (source to relay, relay to destination) = (', num2str(d1(i_d)), ',', num2str(d2(i_d)), ') completed.'])
    disp([' - BER upper bounds, non-CoRe: ', num2str(BER_analytical{i_d}(:, 1)')])
    disp([' - BER upper bounds, Seddik: ', num2str(BER_analytical{i_d}(:, 2)')])
    disp([' - BER upper bounds, QAP: ', num2str(BER_analytical{i_d}(:, 3)')])
    disp([' - BER emperical, non-CoRe: ', num2str(BER_MC{i_d}(:, 1)')])
    disp([' - BER emperical, Seddik: ', num2str(BER_MC{i_d}(:, 2)')])
    disp([' - BER emperical, QAP: ', num2str(BER_MC{i_d}(:, 3)')])
end

BER_analytical = reshape(cell2mat(BER_analytical), M, 3 * n_d);
BER_MC = reshape(cell2mat(BER_MC), M, 3 * n_d);

%% Visualization
% The BER upperbound
cmap = [0, 0, 0 ;0.5, 0, 1; 0, 0, 1; 1, 0, 0];
legend_item = cell(3 * M - 2, 1);
h = figure;
semilogy(d1, (BER_analytical(1, 1 : n_d) + BER_analytical(1, n_d : -1 : 1)) / 2, 'k+-','linewidth', 2), hold on;
legend_item{1} = 'TR1';
for m = 2 : M
    semilogy(d1, (BER_analytical(m, 1 : n_d) + BER_analytical(m, n_d : - 1 : 1 )) / 2, '+-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(d1, (BER_analytical(m, (n_d + 1) : (2 * n_d)) + BER_analytical(m, (2 * n_d) : -1 : (n_d + 1))) / 2, '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(d1, (BER_analytical(m, 2 * n_d + 1 : 3 * n_d) + BER_analytical(m, (3 * n_d) : -1 : 2 * n_d + 1)) / 2, 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 4} = ['NM', num2str(m)];
    legend_item{3 * m - 3} = ['GS', num2str(m)];
    legend_item{3 * m - 2} = ['QAP', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('d_{S-R}'), ylabel('Avg. BER');
xlim([0, 1.5]), ylim([1e-5, 1]);
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_distance_upperbound_', num2str(Q), 'QAM.fig']);

%% The empirical BER
h = figure;
semilogy(d1, (BER_MC(1, 1 : n_d) + BER_MC(1, n_d : -1 : 1)) / 2, 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR1';
for m = 2 : M
    semilogy(d1, (BER_MC(m, 1 : n_d) + BER_MC(m, n_d : -1 : 1)) / 2, '+-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(d1, (BER_MC(m, n_d + 1 : 2 * n_d) + BER_MC(m, (2 * n_d) : -1 : (n_d + 1))) / 2, '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(d1, (BER_MC(m, 2 * n_d + 1 : 3 * n_d) + BER_MC(m, (3 * n_d) : -1 : (2 * n_d + 1))) / 2, 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 4} = ['NM', num2str(m)];
    legend_item{3 * m - 3} = ['GS', num2str(m)];
    legend_item{3 * m - 2} = ['QAP', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('d_{S-R}'), ylabel('Avg. BER');
xlim([0, 1.5]), ylim([1e-5, 1]);
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_distance_MonteCarlo_', num2str(Q), 'QAM.fig']);