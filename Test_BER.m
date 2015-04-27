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
dB_inv_sigma2 = 25; % 1/sigma2 in dB
Pr = 2; % Power at the relay
d1 = 0.5; % Distance between S and R
d2 = 0.5; % Distance between R and D
nu = 3; % Pathloss factor

N = 5000; % Number of monte-carlo run

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

% map = repmat(1 : Q, 2, 1); % The mapping from indices to constellation points (all gray mapping)
map = zeros(2, Q);
map(1, :) = 1 : Q; % Gray mapping
map(2, :) = 1 : Q; % Gray mapping
%map(2, :) = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % I kind of forget where does this mapping come from

%% 3. Not let us test the bit error rate

% Compute the bit error rate using our analytical upper bound
BER_analytical = get_BER_upper_bound(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)

% Compute the bit error rate using Monte-Carlo simulation
BER_MC = get_BER(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, N)
