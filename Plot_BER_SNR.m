clear all;
close all;
clc;

addpath('../functions');

%% 0. Load the data file that contains the test result
load('Test_2015514162030429.mat')

%% 1. Simulation settings
N = 2e7; % Number of monte-carlo run

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
pwr = 1; 

Pr = test_cases(1).param_origin.Pr; % Power at the relay
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
constellation = get_constellation(Nbps, type_mod, pwr);

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