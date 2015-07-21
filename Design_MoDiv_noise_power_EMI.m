% The script to compute and save the 4-D cost matrix

clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% Constellation specification
Nbps = 6; % 4, 5, 6
type_mod = 'QAM';

% Node S, R, D power, channel power and noise power specification
% For now we assume that:
%   S and D are using the same unit power 
%   all the 3 node S, R, D's AWGN noise have the same power sigma2
% We assume that the S and D are separated by distance 1. The variance of
% the relay channel is d ^ -nu where nu is the pathloss factor
% We also assume that the channel and noise are stationary across
% transmissions

dB_inv_sigma2 = [-20 : 2 : 16]; % 1/sigma2 in dB
p_Pr = 0.5; % this portion of the total power of 4 is allocated to the relay. The rest are divided eqaully between the 2 end nodes
d = [0.5, 0.5]; % Distance between S and R, R and D

nu = 3; % Pathloss factor
M = 4; % Number of retransmission

epsilon = 0.01; % Tolerance to control the error of scaling the 2 cost matrices to integer
n_itr = 100000; % Number of iterations for the tabu QAP solver

%% 2. Initialization: generate and save all test cases
test_cases = construct_test_cases(Nbps, type_mod, dB_inv_sigma2, p_Pr, d, nu, M, true);
n_case = length(test_cases);
time_step = regexprep(num2str(clock),'[^\w'']',''); % The time step used to label all saved files as a suffix

%% 3. Start the computation
% Compute the expected number of pairwise error bit based on the cost 
% matrix for the previous transmissions, or initialize this value to be 1 /
% 2 of the hamming distance matrix

for i_case = 1 : n_case
    disp(['Test case ', num2str(i_case), '/', num2str(n_case)]);
    
    % unpack the parameters
    
    E = test_cases(i_case).param_derived.E;
    prod_E_over_Q = E / test_cases(i_case).param_derived.Q; % Initialize the term inside the second log function: the product of E across all retransmissions over Q
    
    pow10_E = get_scale_power10(E, epsilon);
    pow10_prod_E_over_Q = get_scale_power10(prod_E_over_Q, epsilon);
    
    disp([' - Transmission 1: EMI = ', num2str(log2(test_cases(i_case).param_derived.Q) - log2(sum(sum(prod_E_over_Q))))]);
    
    test_cases(i_case).map = zeros(test_cases(i_case).param_origin.M - 1, test_cases(i_case).param_derived.Q);
    for m = 2 : test_cases(i_case).param_origin.M
        
        % Save the two 2-D cost matrices of the QAP-KB problem
        %filename = ['test_case', num2str(i_case), '_', num2str(m) , '_', time_step, '.mat'];
        %save(filename, 'xpcd_PBER', 'E');

        % Put the QAP solver here, currently we use a place holder which
        % returns gray mapping only. Save the resulting map also to the
        % structure
  
        map = solve_QAP_tabu(prod_E_over_Q, E, pow10_prod_E_over_Q, pow10_E, n_itr);
        test_cases(i_case).map(m - 1, :) = map;
        
        % Update the product of E across all retransmissions over Q
        prod_E_over_Q = get_prod_E_over_Q(prod_E_over_Q, E, map);
        pow10_prod_E_over_Q = get_scale_power10(prod_E_over_Q, epsilon);
        % toc;
        disp([' - Transmission ', num2str(m), ': EMI = ', num2str(log2(test_cases(i_case).param_derived.Q) - log2(sum(sum(prod_E_over_Q))))]);
    end
    
    % Save the test results regularly in case of a crash
    save(['Test_', time_step, '.mat'], 'test_cases');
    disp(['Test case ', num2str(i_case), '/', num2str(n_case), ' saved.']);
end
