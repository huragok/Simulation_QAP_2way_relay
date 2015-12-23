clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% 64 QAM, dB_inv_sigma2 = 3 : 2.5 : 7
load('Test_20151222182339679.mat');

% M_to_test = 2; % 3, 4, 5 plot the waterfall curve only for Chase combining M_to_test transmissions
% dB_inv_sigma2 = {[4.5, 5, 5.5, 6, 6.3, 6.5, 6.7, 6.9, 7, 7.1],...
%                  [4.5, 5, 5.5, 6, 6.3, 6.5 ,6.6, 6.7, 6.8, 6.9],...
%                  [4.5, 5, 5.5, 6, 6.3, 6.5, 6.6, 6.7, 6.8, 6.9],...
%                  [4.5, 5, 5.5, 6, 6.3, 6.5, 6.6, 6.7, 6.8, 6.9],...
%                  [4.5, 5, 5.5, 6, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8 ]};

% M_to_test = 3;
% dB_inv_sigma2 = {[1, 1.5, 2, 2.5, 2.8, 3, 3.1, 3.2, 3.3],...
%                  [1, 1.5, 2, 2.5, 2.8, 2.9, 3, 3.1, 3.2],...
%                  [1, 1.5, 2, 2.5, 2.7, 2.8, 2.9, 3, 3.1],...
%                  [1, 1.5, 2, 2.4, 2.6, 2.8, 2.9, 3],...
%                  [1, 1.5, 2, 2.3, 2.5, 2.6, 2.7, 2.8, 2.9]};

% M_to_test = 4;
% dB_inv_sigma2 = {[-1, -0.5, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],...
%                  [-1, -0.5, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],...
%                  [-1, -0.5, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],...
%                  [-1, -0.5, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],...
%                  [-1, -0.5, 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]};

M_to_test = 5;
dB_inv_sigma2 = {[-2.5, -2, -1.5, -1.2, -1, -0.9, -0.8, -0.75],...
                 [-2.5, -2, -1.6, -1.4, -1.2, -1.1, -1, -0.9],...
                 [-2.5, -2, -1.6, -1.4, -1.2, -1.1, -1, -0.9],...
                 [-2.5, -2, -1.6, -1.4, -1.2, -1, -0.9, -0.8, -0.75],...
                 [-2.5, -2, -1.6, -1.4, -1.2, -1, -0.9, -0.8, -0.6]};


Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
d = [test_cases(1).param_origin.d1, test_cases(1).param_origin.d2]; % Distance between S and R, R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

max_frame = 2000;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

%% 2. Initialization
Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

beta_sr = d(1) ^ -nu;
beta_rd = d(2) ^ -nu;

n_case = length(test_cases);
n_sigma2 = zeros(n_case, 1);
map  = cell(n_case, 1);
codedBER = cell(n_case, 1);
for i_case = 1 : n_case
    n_sigma2(i_case) = length(dB_inv_sigma2{i_case});
    map{i_case} = [1 : Q; test_cases(i_case).map];
    codedBER{i_case} = zeros(n_sigma2(i_case), 1);
end

%% 3. Now let us test the coded bit error rate
for i_case = 1 : n_case
    sigma_sqr = 10 .^ (-dB_inv_sigma2{i_case} / 10); % The noise covariance at all nodes
    for i_sigma2 = 1 : n_sigma2(i_case)
        tic
        % Compute the bit error rate using Monte-Carlo simulation
        codedBER{i_case}(i_sigma2) = get_codedBER(constellation, map{i_case}(1 : M_to_test, :), beta_sr, beta_rd, Pr, Pt, Pt, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2), max_frame, iter_max, coding_rate, nldpc, seed);
        toc;
        disp(['Design 1/sigma2 = ', num2str(test_cases(i_case).param_origin.dB_inv_sigma2), 'dB, actual 1/sigma2 = ', num2str(dB_inv_sigma2{i_case}(i_sigma2)), 'dB, coded BER = ', num2str(codedBER{i_case}(i_sigma2))]);
    end
end


%% 4.Visualization
cmap = [0, 0, 0; 0, 0, 1; 1, 0, 0; 0, 1, 0; 1, 1, 0];
h = figure;
for i_case = 1 : n_case
    semilogy(dB_inv_sigma2{i_case}, codedBER{i_case}, '+-', 'Color', cmap(i_case, :), 'linewidth', 2), hold on;
end

save(['waterfall_mismatch_', num2str(M_to_test), 'M_', num2str(Q), 'QAM.mat'], 'dB_inv_sigma2', 'codedBER')
