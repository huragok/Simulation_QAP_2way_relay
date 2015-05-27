clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
Nbps = 4; % 4, 5, 6
type_mod = 'QAM';

dB_inv_sigma2 = 10; % 1/sigma2 in dB
p_Pr = 0.5; % this portion of the total power of 4 is allocated to the relay. The rest are divided eqaully between the 2 end nodes
d = [0.5, 0.5]; % Distance between S and R, R and D

nu = 3; % Pathloss factor
M = 3; % Total number of transmissions

max_frame = 100;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 17;
%% 2. Initialization
n_sigma2 = length(dB_inv_sigma2);

Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d(1) ^ -nu;
beta_rd = d(2) ^ -nu;

g = sqrt(Pr ./ (beta_sr * Pt + beta_rd * Pt + sigma_sqr_r)); % The power normalization factor

map_seddik = get_map_seddik(Q, M);

%% 3. Not let us test the coded bit error rate
BER_coded = zeros(1, n_sigma2);

for i_sigma2 = 1 : n_sigma2
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    BER_coded(i_sigma2) = get_codedBER(constellation, map_seddik, beta_sr, beta_rd, g(i_sigma2), sigma_sqr_d(i_sigma2), sigma_sqr_r(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
  
    toc;
    disp(['Coded BER simulation for 1/sigma2 = ', num2str(dB_inv_sigma2(i_sigma2)), 'dB completed.'])
    disp([' - LDPC coded BER, Seddik: ', num2str(BER_coded(:, i_sigma2)')]);
end
