function E = get_factor_PEP_update_MC(dist_sqr, beta_sr, beta_rd, Pt, Pr, sigma_sqr_d, sigma_sqr_r, N_batch, N_per_batch, seed)
%   E = get_factor_PEP_update_MC(dist_sqr, beta, g, sigma_sqr, sigma_sqr_r, N)
%   Get the factor to be multiplied to PEP between indices p and q for the
%   first (M-1) retransmissions to compute the PEP for the first M 
%   retransmissions. Compute with Monte-Carlo simulation of N run
% _________________________________________________________________________
%	Inputs:
% 		dist_sqr:       Vector, the square norm of the 2 comstellation 
%                       points to which p and q are mapped in the M-th
%                       retransmission
%       beta_sr:        Scalar, the variance of the Rayleigh channel from
%                       source to relay
%       beta_rd:        Scalar, the variance of the Rayleigh channel from
%                       relay to destination
%       sigma_sqr_d:    Scalar, the variance of AWGN noise at the
%                       destination
%       sigma_sqr_r:    Scalar, the variance of AWGN noise at the relay
%       N:              Scalar, size of Monte-Carlo simulation
%	Outputs:
%		E:              Scalar, the factor used to update PEP for the first
%                       (M-1) retransmissions into PEP for the first M
%                       retransmissions between p and q
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/22/2015
% Codename: Dunkirk
% _________________________________________________________________________

[nrow, ncol] = size(dist_sqr);
E = zeros(nrow, ncol);

rng(seed);
for i_batch = 1 : N_batch
    delta1 = abs(sqrt(beta_sr / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    gamma2 = abs(sqrt(beta_rd / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) .^ 2;
    alpha = sqrt(Pr ./ (delta1 * Pt + gamma2 * Pt + sigma_sqr_r)); % The power normalization factor
    sigma_sqr_tilde = sigma_sqr_d + sigma_sqr_r * alpha .^ 2 .* gamma2; % The effective noise power at the destination
    
    tmp = alpha .^ 2 .* gamma2 .* delta1 ./ (sigma_sqr_tilde .^ 2) / 4;
    for irow = 1 : nrow
        for icol = 1 : ncol
            E(irow, icol) = E(irow, icol) + mean(exp(-dist_sqr(irow, icol) * tmp));
        end
    end
end

E = E / N_batch;