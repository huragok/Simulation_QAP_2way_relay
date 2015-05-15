clear all;
close all;
clc;

addpath('./functions');

%% 0. Load the data file that contains the test result
load('Test_2015515132839775.MAT') % 16QAM

%% 1. Simulation settings
N = 1E6; % Number of monte-carlo run

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
dB_inv_sigma2 = test_cases(1).param_origin.dB_inv_sigma2;
d1 = test_cases(1).param_origin.d1;
d2 = test_cases(1).param_origin.d2;
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

n_p_Pr = length(test_cases);
p_Pr = zeros(n_p_Pr, 1);
for i_p_Pr = 1 : n_p_Pr
    p_Pr(i_p_Pr) = test_cases(i_p_Pr).param_origin.p_Pr;
end

%% 2. Initialization
Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 .^ -nu;
beta_rd = d2 .^ -nu;

g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor

map_noncore = get_map_noncore(Q, M);
map_seddik = get_map_seddik(Q, M);
map_QAP = cell(n_p_Pr, 1);
for i_p_Pr = 1 : n_p_Pr
    map_QAP{i_p_Pr} = [1 : Q; test_cases(i_p_Pr).map];
end

%% 3. Not let us test the bit error rate
BER_analytical = cell(n_p_Pr, 1);
BER_MC = cell(n_p_Pr, 1);

for i_p_Pr = 1 : n_p_Pr
    tic

    % Compute the bit error rate using our analytical upper bound
    BER_analytical{i_p_Pr} = zeros(M, 3);
    BER_analytical{i_p_Pr}(:, 1) = get_BER_upper_bound(sqrt(Pt(i_p_Pr)) * constellation, map_noncore, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r);
    BER_analytical{i_p_Pr}(:, 2) = get_BER_upper_bound(sqrt(Pt(i_p_Pr)) * constellation, map_seddik, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r);
    BER_analytical{i_p_Pr}(:, 3) = get_BER_upper_bound(sqrt(Pt(i_p_Pr)) * constellation, map_QAP{i_p_Pr}, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r);
    
    % Compute the bit error rate using Monte-Carlo simulation
    BER_MC{i_p_Pr} = zeros(M, 3);
    %BER_MC{i_p_Pr}(:, 1) = get_BER(sqrt(Pt(i_p_Pr)) * constellation, map_noncore, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r, N);
    %BER_MC{i_p_Pr}(:, 2) = get_BER(sqrt(Pt(i_p_Pr)) * constellation, map_seddik, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r, N);
    %BER_MC{i_p_Pr}(:, 3) = get_BER(sqrt(Pt(i_p_Pr)) * constellation, map_QAP{i_p_Pr}, beta_sr, beta_rd, g(i_p_Pr), sigma_sqr_d, sigma_sqr_r, N);
    
    toc;
    disp(['Power allocated to the relay = ', num2str(p_Pr(i_p_Pr)), ' completed.'])
    disp([' - BER upper bounds, non-CoRe: ', num2str(BER_analytical{i_p_Pr}(:, 1)')])
    disp([' - BER upper bounds, Seddik: ', num2str(BER_analytical{i_p_Pr}(:, 2)')])
    disp([' - BER upper bounds, QAP: ', num2str(BER_analytical{i_p_Pr}(:, 3)')])
    disp([' - BER emperical, non-CoRe: ', num2str(BER_MC{i_p_Pr}(:, 1)')])
    disp([' - BER emperical, Seddik: ', num2str(BER_MC{i_p_Pr}(:, 2)')])
    disp([' - BER emperical, QAP: ', num2str(BER_MC{i_p_Pr}(:, 3)')])
end

BER_analytical = reshape(cell2mat(BER_analytical), M, 3 * n_p_Pr);
BER_MC = reshape(cell2mat(BER_MC), M, 3 * n_p_Pr);

%% Visualization
% The BER upperbound
cmap = colormap(hsv(M));
legend_item = cell(3 * M, 1);
h = figure;
for m = 1 : M
    semilogy(p_Pr, BER_analytical(m, 1 : n_p_Pr), '-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(p_Pr, BER_analytical(m, n_p_Pr + 1 : 2 * n_p_Pr), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(p_Pr, BER_analytical(m, 2 * n_p_Pr + 1 : 3 * n_p_Pr), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 2} = ['Non-CoRe, M = ', num2str(m)];
    legend_item{3 * m - 1} = ['Seddik, M = ', num2str(m)];
    legend_item{3 * m} = ['MoDiv, M = ', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('p_R'), ylabel('BER');
legend(legend_item, 'Location', 'eastoutside');
saveas(h, ['BER_power_allocation_upperbound_', num2str(Q), 'QAM.fig']);

%% The empirical BER
h = figure;
for m = 1 : M
    semilogy(p_Pr, BER_MC(m, 1 : n_p_Pr), '-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(p_Pr, BER_MC(m, n_p_Pr + 1 : 2 * n_p_Pr), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(p_Pr, BER_MC(m, 2 * n_p_Pr + 1 : 3 * n_p_Pr), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 2} = ['Non-CoRe, M = ', num2str(m)];
    legend_item{3 * m - 1} = ['Seddik, M = ', num2str(m)];
    legend_item{3 * m} = ['MoDiv, M = ', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('p_R'), ylabel('BER');
legend(legend_item, 'Location', 'eastoutside');
saveas(h, ['BER_power_allocation_MonteCarlo_', num2str(Q), 'QAM.fig']);