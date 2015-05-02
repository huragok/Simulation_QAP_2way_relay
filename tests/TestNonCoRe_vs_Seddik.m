% Script to test the compare the BER upperbound and the Monte-Carlo 
% simulated BER of different remapping schemes. 
% Based on script "Test_BER.m".

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
dB_inv_sigma2 = [16 : 2 : 30]; % 1/sigma2 in dB
Pr = 2; % Power at the relay
d1 = 0.5; % Distance between S and R
d2 = 0.5; % Distance between R and D
nu = 3; % Pathloss factor

N = 2e7; % Number of monte-carlo run

M = 4; % Total number of transmissions

sedik = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1;
%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, pwr);

n_sigma2 = length(dB_inv_sigma2);
sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr ./ (beta_sr + beta_rd + sigma_sqr_r)); % The power normalization factor

map_noncore = repmat(1 : Q, M, 1);
map_seddik = zeros(M, Q);
map_seddik(1, :) = 1 : Q; % Gray mapping
map_seddik(2 : M, :) = repmat(sedik, M - 1, 1);

%% 3. Not let us test the bit error rate

BER_analytical = cell(n_sigma2, 1);
BER_MC = cell(n_sigma2, 1);
for i_sigma2 = 1 : n_sigma2
    tic
    % Compute the bit error rate using our analytical upper bound
    BER_analytical{i_sigma2} = zeros(M, 2);
    BER_analytical{i_sigma2}(:, 1) = get_BER_upper_bound(constellation, map_noncore, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    BER_analytical{i_sigma2}(:, 2) = get_BER_upper_bound(constellation, map_seddik, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2));
    
    % Compute the bit error rate using Monte-Carlo simulation
    BER_MC{i_sigma2} = zeros(M, 2);
    BER_MC{i_sigma2}(:, 1) = get_BER(constellation, map_noncore, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N);
    BER_MC{i_sigma2}(:, 2) = get_BER(constellation, map_seddik, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N);
    toc;
    disp(['BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - BER upper bounds, non-CoRe: ', num2str(BER_analytical{i_sigma2}(:, 1)')])
    disp([' - BER emperical, non-CoRe: ', num2str(BER_MC{i_sigma2}(:, 1)')])
    disp([' - BER upper bounds, Seddik: ', num2str(BER_analytical{i_sigma2}(:, 2)')])
    disp([' - BER emperical, Seddik: ', num2str(BER_MC{i_sigma2}(:, 2)')])
end

BER_analytical = reshape(cell2mat(BER_analytical), M, 2 * n_sigma2);
BER_MC = reshape(cell2mat(BER_MC), M, 2 * n_sigma2);

%% Visualization
cmap = colormap(hsv(M));
legend_item = cell(4 * M, 1);
h = figure;
for m = 1 : M
    semilogy(dB_inv_sigma2, BER_analytical(m, 1 : n_sigma2), 'o-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_MC(m, 1 : n_sigma2), 'o--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_analytical(m, n_sigma2 + 1 :  2 * n_sigma2), '+-', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    semilogy(dB_inv_sigma2, BER_MC(m, n_sigma2 + 1 :  2 * n_sigma2), '+--', 'Color', cmap(m, :), 'linewidth', 2), hold on;
    legend_item{4 * m - 3} = ['Non-CoRe bound, M = ', num2str(m)];
    legend_item{4 * m - 2} = ['Non-CoRe MC, M = ', num2str(m)];
    legend_item{4 * m - 1} = ['Seddik bound, M = ', num2str(m)];
    legend_item{4 * m} = ['Seddik MC, M = ', num2str(m)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
legend(legend_item);
saveas(h, 'Test_NonCoRe_vs_Seddik.fig');

