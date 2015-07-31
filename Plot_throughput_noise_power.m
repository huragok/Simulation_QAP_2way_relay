clear all;
close all;
clc;

addpath('./functions');

load('Test_201573023658.mat')
M_to_test = 6;
dB_inv_sigma2_noncore = [12];
dB_inv_sigma2_core = [12];
dB_inv_sigma2_QAP = [12];

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
