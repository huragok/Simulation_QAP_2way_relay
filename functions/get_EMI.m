function EMI = get_EMI(constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed)
%   EMI = get_EMI(constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed)
%   Get the EMI by actually running Monte-Carlo simulation to compute the
%   expectation for EMI numerically
% _________________________________________________________________________
%	Inputs:
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta_sr:        Scalar, the variance of the Rayleigh channel from
%                       source to relay
%       beta_rd:        Scalar, the variance of the Rayleigh channel from
%                       relay to destination
%       Pr:             Scalar, the average power constraint at the relay
%       P1:             Scalar, the average power constraint at the
%                       source
%       P2:             Scalar, the average power constraint at the
%                       destination
%       sigma_sqr_d:    Scalar, the variance of AWGN noise at the
%                       destination
%       sigma_sqr_r:    Scalar, the variance of AWGN noise at the relay
%       N_per_batch:    Scalar, number of Monte-Carlo run per batch (size 
%                       of vectorization)
%       N_batch:        Scalar, number of batches (for-loop size)
%       seed:           Scalar, seed for the random number generator
%	Outputs:
%		EMI:			M-by-1 vector, the EMI after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 07/20/2015
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map);
symbols_mapped = constellation(map); % The mapped symbols at all transmissions

EMI = zeros(M, N_batch);

rng(seed);
for i_batch = 1 : N_batch

    exp_nnc = ones(Q, N_per_batch, Q);
    for m = 1 : M
        h_sr = sqrt(beta_sr / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
        h_rd = sqrt(beta_rd / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));
        g = sqrt(Pr ./ (abs(h_sr) .^ 2 * P1 + abs(h_rd) .^ 2 * P2 + sigma_sqr_r)); % The power normalization factor
        sigma = sqrt(sigma_sqr_d + sigma_sqr_r * g .^ 2 .* h_rd .^ 2);
        n = sqrt(1 / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));
        for p = 1 : Q
            exp_nnc(:, :, p) = exp_nnc(:, :, p) .* exp(repmat(abs(n) .^ 2, Q, 1) - abs(repmat(n, Q, 1) - (symbols_mapped(m, p) - symbols_mapped(m, :)).' * (g .* h_rd .* h_sr ./ sigma)) .^ 2);
        end

        EMI(m, i_batch) = log2(Q) - 1 / Q * mean(sum(log2(sum(exp_nnc, 1)), 3)); % EMI for the m-th transmission
    end
end

EMI = mean(EMI, 2);

