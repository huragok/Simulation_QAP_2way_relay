% Script to test the functionality of get_factor_PEP_update function and
% compare its result with the true factor in the form of E[exp()] using
% Monte-Carlo simulation

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

N = 5000; % Number of monte-carlo run

%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, pwr);

sigma_sqr = 10 ^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr / (beta_sr + beta_rd + sigma_sqr_r)); % The power normalization factor

%% 3. Test some  E
dist_sqr = (0 : 0.1 : 1);
E = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
