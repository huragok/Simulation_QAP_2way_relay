clear all;
close all;
clc;

addpath('./functions');

% 64QAM, retransmission 2 and 3 6 : 2 : 20
% load('Test_201512221170973.mat'); % E evaluated with Monte-Carlo simulation
% test_cases_MC = test_cases;
% load('Test_2015916152041349.mat'); % E evaluated with Proposition 1


% 64QAM, retransmission 4 and 5 1 : 1.5 : 13
% load('Test_20151222115413653.mat'); % E evaluated with Monte-Carlo simulation
% test_cases_MC = test_cases;
% load('Test_201591618110915.mat')

% 16QAM, retransmission 2 and 3 2 : 2 : 16
% load('Test_2015122214632161.mat'); % E evaluated with Monte-Carlo simulation
% test_cases_MC = test_cases;
% load('Test_2015916132547305.mat'); % E evaluated with Proposition 1

% 16QAM, retransmission 4 and 5 -2 : 1.5 : 10
load('Test_2015122214168195.mat'); % E evaluated with Monte-Carlo simulation
test_cases_MC = test_cases;
load('Test_2015916143412484.mat'); % E evaluated with Proposition 1

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

map_QAP = cell(n_sigma2, 1);
map_QAP_MC = cell(n_sigma2, 1);
for i_sigma2 = 1 : n_sigma2
    map_QAP{i_sigma2} = [1 : Q; test_cases(i_sigma2).map];
    map_QAP_MC{i_sigma2} = [1 : Q; test_cases_MC(i_sigma2).map];
end

%% 3. Now let us test the bit error rate
BER = zeros(M, n_sigma2); % Simulated BER with E evaluated analytically
BER_MC = zeros(M, n_sigma2);  % Simulated BER with E evaluated using Monte-Carlo simulation

for i_sigma2 = 1 : n_sigma2
    tic
    
    BER(:, i_sigma2) = get_BER(constellation, map_QAP{i_sigma2}, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    BER_MC(:, i_sigma2) = get_BER(constellation, map_QAP_MC{i_sigma2}, beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), N_per_batch, N_batch, seed);
    
    toc;
    disp(['BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - BER using analytical E: ', num2str(BER(:, i_sigma2)')])
    disp([' - BER using numerical E: ', num2str(BER_MC(:, i_sigma2)')])
end

save(['BER_noise_power_MC_', num2str(Q), 'QAM.mat'], 'dB_inv_sigma2', 'BER', 'BER_MC');

%% Visualization
%% The empirical BER
cmap = [0, 0, 0; 0, 0, 1; 1, 0, 0; 0, 1, 0; 1, 1, 0];
legend_item = cell(2 * M - 1, 1);
h = figure;
semilogy(dB_inv_sigma2, BER(1, 1 : n_sigma2), 'k+-', 'linewidth', 2), hold on;
legend_item{1} = 'TR0';
for m = 2 : M
    semilogy(dB_inv_sigma2, BER(m, :), '+-', 'Color', cmap(m, :), 'linewidth', 2);
    semilogy(dB_inv_sigma2, BER_MC(m, :), '^--', 'Color', cmap(m, :), 'linewidth', 2), hold on;

    legend_item{2 * m - 2} = ['Approx', num2str(m-1)];
    legend_item{2 * m - 1} = ['Num', num2str(m-1)];
end
grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('BER');
ylim([1e-6, 10]), xlim([2, 28])
legend(legend_item, 'Location', 'northeast');
saveas(h, ['BER_noise_power_MC_MonteCarlo_', num2str(Q), 'QAM.fig']);