clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
load('Test_2015526214647641.MAT') % 16QAM
M_to_test = 3; % 3, 4, plot the waterfall curve only for Chase combining M_to_test transmissions
dB_inv_sigma2_noncore = [-1, -0.5]; % 1/sigma2 in dB
dB_inv_sigma2_seddik = [-1, -0.5];
dB_inv_sigma2_QAP = [-1, -0.5];


Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
d = [test_cases(1).param_origin.d1, test_cases(1).param_origin.d2]; % Distance between S and R, R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

max_frame = 100;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

%% 2. Initialization
n_sigma2_noncore = length(dB_inv_sigma2_noncore);
n_sigma2_seddik = length(dB_inv_sigma2_seddik);
n_sigma2_QAP = length(dB_inv_sigma2_QAP);

Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

beta_sr = d(1) ^ -nu;
beta_rd = d(2) ^ -nu;

map_noncore = get_map_noncore(Q, M);
map_seddik = get_map_seddik(Q, M);
map_QAP  = [1 : Q; test_cases(1).map];

%% 3. Not let us test the coded bit error rate
codedBER_noncore = zeros(1, n_sigma2_noncore); % coded BER for noncore scheme
codedBER_seddik = zeros(1, n_sigma2_seddik);
codedBER_QAP = zeros(1, n_sigma2_QAP);

% The non-CoRe waterfall curve
sigma_sqr = 10 .^ (-dB_inv_sigma2_noncore / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor
for i_sigma2 = 1 : n_sigma2_noncore
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_noncore(i_sigma2) = get_codedBER(constellation, map_noncore(1 : M_to_test, :), beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Non-CoRe, 1/sigma2 = ', num2str(dB_inv_sigma2_noncore(i_sigma2)), 'dB, coded BER = ', num2str(codedBER_noncore(i_sigma2))]);
end

% The Seddik waterfall curve
sigma_sqr = 10 .^ (-dB_inv_sigma2_seddik / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor
for i_sigma2 = 1 : n_sigma2_seddik
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_seddik(i_sigma2) = get_codedBER(constellation, map_seddik(1 : M_to_test, :), beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Seddik, 1/sigma2 = ', num2str(dB_inv_sigma2_seddik(i_sigma2)), 'dB, coded BER = ', num2str(codedBER_seddik(i_sigma2))]);
end

% The QAP waterfall curve
sigma_sqr = 10 .^ (-dB_inv_sigma2_QAP / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor
for i_sigma2 = 1 : n_sigma2_QAP
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_QAP(i_sigma2) = get_codedBER(constellation, map_QAP(1 : M_to_test, :), beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['QAP, 1/sigma2 = ', num2str(dB_inv_sigma2_QAP(i_sigma2)), 'dB, coded BER = ', num2str(codedBER_QAP(i_sigma2))]);
end

%% 4. Visualization
% The BER upperbound
h = figure;
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'k+-', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_seddik, codedBER_seddik, 'b^--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'ro-.', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
legend({'non-Core', 'Seddik', 'MoDiv'}, 'Location', 'eastoutside');
saveas(h, ['waterfall_', num2str(M_to_test), 'M_', num2str(Q), 'QAM.fig']);


