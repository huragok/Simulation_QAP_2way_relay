% The script to compute and save the 4-D cost matrix

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


%% 2. Initialization
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, pwr);

sigma_sqr = 10 ^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes
sigma_sqr_d = sigma_sqr;
sigma_sqr_r = sigma_sqr;

beta_sr = d1 ^ -nu;
beta_rd = d2 ^ -nu;

g = sqrt(Pr / (beta_sr + beta_rd + sigma_sqr_r)); % The power normalization factor
%% 3. Start the computation
% Compute the expected number of pairwise error bit based on the cost 
% matrix for the previous transmissions, or initialize this value to be 1 /
% 2 of the hamming distance matrix

xpcd_num_PE_bits = get_hamming_dist(Nbps) / 2; % Initialization

%c_last = fscanf('test_last.data', '%f');
%map = 1 : Q; % Gray mapping
%xpcd_num_PE_bits = get_xpcd_num_PE_bits(c, map); % Update 

% Compute the distances betweem each pair of constellation points
dist_sqr = abs(repmat(constellation, 1, Q) - repmat(constellation.', Q, 1)) .^ 2;
dist_sqr = reshape(dist_sqr, Q ^ 2, 1);

E = zeros(Q ^ 2, 1); % The variable to store all Eik
% Main bottleneck, use par for to speed up
for idx_ik = 1 : Q ^ 2
    k = mod(idx_ik - 1, Q) + 1;
    i = floor((idx_ik - 1) / Q) + 1;
    if i ~= k
        E(idx_ik) = get_factor_PEP_update(dist_sqr(idx_ik), beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r);
    end
end

% Compute and save c
c = zeros(1, Q ^ 4);
for idx = 1 : Q ^ 4
    piqk = idx2piqk(idx, Q);
    c(idx) = E((piqk(2) - 1) * Q + piqk(4)) * xpcd_num_PE_bits(piqk(1), piqk(3));
end

fileID = fopen('test.data', 'w+');
fprintf(fileID, '  %18.16e', c);
fclose(fileID);