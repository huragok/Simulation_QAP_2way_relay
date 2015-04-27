% Script to test the BER upperbound based on the QAP formulation and the
% compare it with the Monte-Carlo simulated BER,

clear all;
close all;
clc;

addpath('./functions');

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
dB_inv_sigma2 = 0; % 1/sigma2 in dB
Pr = 2; % Power at the relay
d1 = 0.5; % Distance between S and R
d2 = 0.5; % Distance between R and D
nu = 3; % Pathloss factor

N = 50000; % Number of monte-carlo run

M = 2; % Total number of transmissions
%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, pwr);

sigma_sqr = 10 ^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr / (beta_sr + beta_rd + sigma_sqr_r)); % The power normalization factor

h_sr = sqrt(beta_sr / 2) * (randn(N, M) + 1i * randn(N, M)); % Generate the Rayleigh channel 
h_rd = sqrt(beta_rd / 2) * (randn(N, M) + 1i * randn(N, M));

n_r = sqrt(sigma_sqr_r / 2) * (randn(N, M) + 1i * randn(N, M)); % Generate the random noise
n_d = sqrt(sigma_sqr_d / 2) * (randn(N, M) + 1i * randn(N, M)); % Generate the random noise

% map = repmat(1 : Q, 2, 1); % The mapping from indices to constellation points (all gray mapping)
map = zeros(2, Q);
map(1, :) = 1 : Q; % Gray mapping
map(2, :) = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % I kind of forget where does this mapping come from

%% 3. Not let us test the bit error rate

% Compute the bit error rate using our analytical upper bound
xpcd_num_PE_bits = get_hamming_dist(Nbps) / 2 / Q; % Initialization
for m = 1 : M
    dist_sqr = abs(repmat(constellation(map(m, :)), 1, Q) - repmat(constellation(map(m, :)).', Q, 1)) .^ 2; % A Q-by-Q matrix containing the distance square mesurements
    E = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r); % Get this thing fully vectorized
    xpcd_num_PE_bits = xpcd_num_PE_bits .* E; % Update the expected number of pairwise error beat
end
BER_analytical = sum(sum(xpcd_num_PE_bits));

% Compute the bit error rate using Monte-Carlo simulation
