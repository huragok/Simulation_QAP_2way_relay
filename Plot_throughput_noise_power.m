clear all;
close all;
clc;

addpath('./functions');

% load('Test_201573023658.mat')
% M_to_test = 6;
% dB_inv_sigma2_noncore = [3.5, 4, 4.5, 4.8, 5, 5.5, 6, 6.1, 6.2, 6.5, 7, 7.4, 7.5, 8, 8.5, 9, 9.5, 9.8, 11, 12.5, 12.7, 13.5, 14, 14.2, 14.4];
% dB_inv_sigma2_core = [-1, -0.5, -0.3, 0, 0.4, 0.5, 1, 1.5, 2, 2.2, 3, 3.5, 4, 4.5, 4.8, 5, 6.5, 7, 7.5, 8, 8.3, 11, 12.5, 12.7, 13.5, 14, 14.2, 14.4];
% dB_inv_sigma2_QAP = [-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.2, 2.5, 2.6, 3, 5, 5.5, 6, 6.5, 6.7, 7, 10, 12.5, 12.7, 13.5, 14, 14.2, 14.4];

load('Test_2015914221356302.mat')
M_to_test = 4;
%
dB_inv_sigma2_noncore = [0.2, 0.6, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.7, 6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6];
dB_inv_sigma2_core = [-2, -1.8, -1.6, -0.9, -0.4, 0, 0.3, 0.5, 2.5, 2.8, 3, 3.3, 3.7, 3.9, 6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6];
dB_inv_sigma2_QAP = [-2.4, -2.2, -1.8, -1.6, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 2, 2.2, 2.5, 2.8, 3, 3.5, 5.6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6];

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
d = [test_cases(1).param_origin.d1, test_cases(1).param_origin.d2]; % Distance between S and R, R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

max_frame = 200;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

%% 2. Initialization
n_sigma2_noncore = length(dB_inv_sigma2_noncore);
n_sigma2_core = length(dB_inv_sigma2_core);
n_sigma2_QAP = length(dB_inv_sigma2_QAP);

Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

beta_sr = d(1) ^ -nu;
beta_rd = d(2) ^ -nu;

map_noncore = get_map_noncore(Q, M);
map_core = get_map_core(Q, M);
map_QAP  = [1 : Q; test_cases(1).map];

%% 3. Not let us test the throughput
throughput_noncore = zeros(1, n_sigma2_noncore); % coded BER for noncore scheme
throughput_core = zeros(1, n_sigma2_core);
throughput_QAP = zeros(1, n_sigma2_QAP);

% For non-core
sigma_sqr = 10 .^ (-dB_inv_sigma2_noncore / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
for i_sigma2 = 1 : n_sigma2_noncore
    tic
    % Compute the throughput using Monte-Carlo simulation
    throughput_noncore(i_sigma2) = get_throughput(constellation, map_noncore(1 : M_to_test, :), beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Non-CoRe, 1/sigma2 = ', num2str(dB_inv_sigma2_noncore(i_sigma2)), 'dB, throughput = ', num2str(throughput_noncore(i_sigma2))]);
end

% For core
sigma_sqr = 10 .^ (-dB_inv_sigma2_core / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
for i_sigma2 = 1 : n_sigma2_core
    tic
    % Compute the throughput using Monte-Carlo simulation
    throughput_core(i_sigma2) = get_throughput(constellation, map_core(1 : M_to_test, :), beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['CoRe, 1/sigma2 = ', num2str(dB_inv_sigma2_core(i_sigma2)), 'dB, throughput = ', num2str(throughput_core(i_sigma2))]);
end

% For QAP
sigma_sqr = 10 .^ (-dB_inv_sigma2_QAP / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;
for i_sigma2 = 1 : n_sigma2_QAP
    tic
    % Compute the throughput using Monte-Carlo simulation
    throughput_QAP(i_sigma2) = get_throughput(constellation, map_QAP(1 : M_to_test, :), beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['QAP, 1/sigma2 = ', num2str(dB_inv_sigma2_QAP(i_sigma2)), 'dB, throughput = ', num2str(throughput_QAP(i_sigma2))]);
end


%% 4. Visualization
% The throughput
h = figure;
plot(dB_inv_sigma2_noncore, throughput_noncore, 'k+-', 'linewidth', 2), hold on;
plot(dB_inv_sigma2_core, throughput_core, 'b^--', 'linewidth', 2);
plot(dB_inv_sigma2_QAP, throughput_QAP, 'ro-.', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Throughput (info bits/symbol)');
legend({'non-Core', 'CoRe', 'MoDiv'}, 'Location', 'northeast');
saveas(h, ['throughput_', num2str(M_to_test), 'M_', num2str(Q), 'QAM.fig']);

save(['throughput_', num2str(M_to_test), 'M_', num2str(Q), 'QAM.mat'], 'dB_inv_sigma2_noncore', 'dB_inv_sigma2_core', 'dB_inv_sigma2_QAP', 'throughput_noncore', 'throughput_core', 'throughput_QAP')
