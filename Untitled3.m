clear all;
close all;
clc;

addpath('./functions');


%% 1. Simulation settings
Nbps = 5;
type_mod = 'QAM';

p_Pr = 0.5; % Power at the relay
d1 = 0.5; % Distance between S and R
d2 = 0.5; % Distance between R and D
nu = 3; % Pathloss factor
M = 2; % Total number of transmissions


dB_inv_sigma2 = 16 : 2 : 32;

%% 2. Initialization
n_sigma2 = length(dB_inv_sigma2);
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
map_seddik1 = get_map_seddik(Q, M);
map_seddik2 = get_map_seddik2(Q, M);

%% 3. Not let us test the bit error rate
BER_analytical = cell(n_sigma2, 1);

for i_sigma2 = 1 : n_sigma2
    tic

    % Compute the bit error rate using our analytical upper bound
    BER_analytical{i_sigma2} = zeros(M, 3);
    BER_analytical{i_sigma2}(:, 1) = get_BER_upper_bound(constellation, map_noncore, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    BER_analytical{i_sigma2}(:, 2) = get_BER_upper_bound(constellation, map_seddik1, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    BER_analytical{i_sigma2}(:, 3) = get_BER_upper_bound(constellation, map_seddik2, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    
    
    toc;
    disp(['BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - BER upper bounds, non-CoRe: ', num2str(BER_analytical{i_sigma2}(:, 1)')])
    disp([' - BER upper bounds, Seddik: ', num2str(BER_analytical{i_sigma2}(:, 2)')])
    disp([' - BER upper bounds, QAP: ', num2str(BER_analytical{i_sigma2}(:, 3)')])
end

BER_analytical = reshape(cell2mat(BER_analytical), M, 3 * n_sigma2);

%% Visualization
% The BER upperbound
cmap = colormap(hsv(M));
legend_item = cell(3 * M, 1);
h = figure;
for m = 1 : M
    semilogy(dB_inv_sigma2, BER_analytical(m, 1 : n_sigma2), '-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical(m, n_sigma2 + 1 : 2 * n_sigma2), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical(m, 2 * n_sigma2 + 1 : 3 * n_sigma2), 'o-.', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{3 * m - 2} = ['Non-CoRe, M = ', num2str(m)];
    legend_item{3 * m - 1} = ['Seddik1, M = ', num2str(m)];
    legend_item{3 * m} = ['Seddik2, M = ', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
legend(legend_item, 'Location', 'eastoutside');
saveas(h, ['BER_noise_power_upperbound_', num2str(Q), 'QAM.fig']);