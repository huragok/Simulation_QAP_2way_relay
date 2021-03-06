function test_cases = construct_test_cases(Nbps, type_mod, dB_inv_sigma2, p_Pr, d, nu, M, verbose)
%   test_cases = construct_test_cases(Nbps, type_mod, dB_inv_sigma2, p_Pr, d, nu, M)
%   Construct a structure array for the test cases for the QAP simulation,
%   where the parameter setting is a cartesian product of the input
%   arguments.
% _________________________________________________________________________
%	Inputs:
%       Nbps:           Scalar, number of bits per symbol
%       type_mod:       String, either 'QAM' or 'PSK'
%       dB_inv_sigma2:  n_sigma2-by-1 vector, all 1/sigma2 in dB
%       p_Pr:           n_p_Pr-by-1 vector, the proportion of power
%                       allocated to the relay, assuming the total power 
%                       (2 end nodes and the relay) is 4
%       d:              n_d-by-2 array, the distance between S-R and R-D
%       nu:             Scalar, the pathloss factor
%       M:              Scalar, the total number of retransmissions
%	Outputs:
%       test_cases:     (n_d*n_Pr*n_sigma2)-by-1 structure, each test case
%                       contains all parameters that can fully specify this
%                       test.
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/27/2015
% Codename: Dunkirk
% _________________________________________________________________________

n_sigma2 = length(dB_inv_sigma2);
n_p_Pr = length(p_Pr);
n_d = size(d, 1);

% Derive the common parameters across different cases
Q = 2 ^ Nbps;
constellation = get_constellation(Nbps, type_mod, 1);

sigma_sqr = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes are the same n_sigma2-by-1 vector

Pr = 4 * p_Pr; % Power at the relay
Pt = 2 * (1 - p_Pr); % Power at each of the 2 transmitter

beta_sr = d(:, 1) .^ -nu; % n_d-by-1 vector
beta_rd = d(:, 2) .^ -nu;

test_cases = struct();
i_case = 1;
n_case = n_d * n_p_Pr * n_sigma2;
for i_d = 1 : n_d
    for i_p_Pr = 1 : n_p_Pr
        % Compute the distances betweem each pair of constellation points. We
        % assume that the simulation settings are stationary across transmissions
        dist_sqr = Pt(i_p_Pr) * abs(repmat(constellation, 1, Q) - repmat(constellation.', Q, 1)) .^ 2;

        for i_sigma2 = 1 : n_sigma2
            if verbose
                tic;
            end
            
            % The original parameters
            test_cases(i_case).param_origin.Nbps = Nbps;
            test_cases(i_case).param_origin.type_mod = type_mod;
            test_cases(i_case).param_origin.dB_inv_sigma2 = dB_inv_sigma2(i_sigma2);
            test_cases(i_case).param_origin.p_Pr = p_Pr(i_p_Pr);
            test_cases(i_case).param_origin.d1 = d(i_d, 1);
            test_cases(i_case).param_origin.d2 = d(i_d, 2);
            test_cases(i_case).param_origin.nu = nu;
            test_cases(i_case).param_origin.M = M;
            
            % The derived parameters
            test_cases(i_case).param_derived.Q = Q;
            test_cases(i_case).param_derived.constellation = sqrt(Pt(i_p_Pr)) * constellation;
            test_cases(i_case).param_derived.sigma_sqr_d = sigma_sqr(i_sigma2);
            test_cases(i_case).param_derived.sigma_sqr_r = sigma_sqr(i_sigma2);

            test_cases(i_case).param_derived.beta_sr = beta_sr(i_d);
            test_cases(i_case).param_derived.beta_rd = beta_rd(i_d);

            test_cases(i_case).param_derived.g = sqrt(Pr(i_p_Pr) / (beta_sr(i_d) * Pt(i_p_Pr) + beta_rd(i_d) * Pt(i_p_Pr) + sigma_sqr(i_sigma2))); % The power normalization factor

            % Compute the updating matrix. We assume that the simulation 
            % settings are stationary across transmissions. Saved as a 
            % Q-by-Q matrix
            test_cases(i_case).param_derived.E = get_factor_PEP_update(dist_sqr, beta_sr(i_d), beta_rd(i_d), test_cases(i_case).param_derived.g, sigma_sqr(i_sigma2), sigma_sqr(i_sigma2)); % Get this thing fully vectorized
            
            if verbose
                toc;
                disp(['Test case ', num2str(i_case), '/', num2str(n_case), ' completed.'])
            end
            i_case = i_case + 1;
        end
    end
end
