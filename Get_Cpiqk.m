% The script to compute and save the 4-D cost matrix

clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% Constellation specification
Nbps = 4;
type_mod = 'QAM';

% Node S, R, D power, channel power and noise power specification
% For now we assume that:
%   S and D are using the same unit power 
%   all the 3 node S, R, D's AWGN noise have the same power sigma2
% We assume that the S and D are separated by distance 1. The variance of
% the relay channel is d ^ -nu where nu is the pathloss factor
% We also assume that the channel and noise are stationary across
% transmissions

dB_inv_sigma2 = 16 : 2 : 22; % 1/sigma2 in dB
Pr = 2; % Power at the relay
d = [0.5, 0.5]; % Distance between S and R, R and D

nu = 3; % Pathloss factor
M = 4; % Number of retransmission
%% 2. Initialization: generate and save all test cases
test_cases = construct_test_cases(Nbps, type_mod, dB_inv_sigma2, Pr, d, nu, M, true);
n_case = length(test_cases);
time_step = regexprep(num2str(clock),'[^\w'']',''); % The time step used to label all saved files as a suffix

%% 3. Start the computation
% Compute the expected number of pairwise error bit based on the cost 
% matrix for the previous transmissions, or initialize this value to be 1 /
% 2 of the hamming distance matrix

for i_case = 1 : n_case
    disp(['Test case ', num2str(i_case), '/', num2str(n_case)]);
    
    % unpack the parameters
    
    xpcd_PBER = get_hamming_dist(test_cases(i_case).param_origin.Nbps) / 2 / test_cases(i_case).param_derived.Q / test_cases(i_case).param_origin.Nbps; % Initialize the expected pairwise BER before any transmission
    xpcd_PBER = xpcd_PBER .* test_cases(i_case).param_derived.E; % The expected pairwise BER after the first transmission (Gray mapping)
    disp([' - Transmission 1: BER = ', num2str(sum(sum(xpcd_PBER)))]);
    
    test_cases(i_case).map = zeros(test_cases(i_case).param_origin.M - 1, test_cases(i_case).param_derived.Q);
    for m = 2 : test_cases(i_case).param_origin.M
        % tic;
        % Compute and save the cost matrix for the m-th retransmission
        c = zeros(1, test_cases(i_case).param_derived.Q ^ 4);
        for idx = 1 : test_cases(i_case).param_derived.Q ^ 4
            piqk = idx2piqk(idx, test_cases(i_case).param_derived.Q);
            c(idx) = test_cases(i_case).param_derived.E(piqk(2), piqk(4)) * xpcd_PBER(piqk(1), piqk(3));
        end
        
        % Write each cost matrix to a file
        filename = ['test_case', num2str(i_case), '_', num2str(m) , '_', time_step, '.data'];
        fileID = fopen(filename, 'w+');
        fprintf(fileID, '  %18.16e', c);
        fclose(fileID);

        % Put the QAP solver here, currently we use a place holder which
        % returns gray mapping only. Save the resulting map also to the
        % structure
        map = solve_QAP(c);
        test_cases(i_case).map(m - 1, :) = map;
        
        % Update the expected PBER
        xpcd_PBER = get_xpcd_PBER(c, map);
        % toc;
        disp([' - Transmission ', num2str(m), ': BER = ', num2str(sum(sum(xpcd_PBER)))]);
    end
    
    % Save the test results regularly in case of a crash
    save(['Test_', time_step, '.mat'], 'test_cases');
    disp(['Test case ', num2str(i_case), '/', num2str(n_case), ' saved.']);
end





