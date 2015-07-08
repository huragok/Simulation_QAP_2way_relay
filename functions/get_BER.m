function BER = get_BER(constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed)
%   BER = get_BER(constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, N_per_batch, N_batch, seed)
%   Get the BER by actually running Monte-Carlo simulation on the ML
%   demodulators
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
%		BER:			M-by-1 vector, the BER after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/27/2015
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map);
Nbps = round(log2(Q)); % Number of bit per symbol
symbols_mapped = constellation(map); % The mapped symbols at all transmissions

BER = zeros(M, N_batch);
B = get_hamming_dist(Nbps); % The hamming distance matrix

rng(seed);
for i_batch = 1 : N_batch
    p = randi(Q, 1, N_per_batch); % Generate the random transmitted index
    d = zeros(Q, N_per_batch); % The distance measurement used for the ML demodulator

    for m = 1 : M
        h_sr = sqrt(beta_sr / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
        h_rd = sqrt(beta_rd / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));
        g = sqrt(Pr ./ (abs(h_sr) .^ 2 * P1 + abs(h_rd) .^ 2 * P2 + sigma_sqr_r)); % The power normalization factor
        symbol = symbols_mapped(m, p); % Generate the random symbol

        y = g .* h_rd .* (h_sr .* symbol + sqrt(sigma_sqr_r / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch))) + sqrt(sigma_sqr_d / 2) * (randn(1, N_per_batch) + 1i * randn(1, N_per_batch));

        d = d + abs(repmat(y, Q, 1) - symbols_mapped(m, :).' * (g .* h_rd .* h_sr)) .^ 2 ./ repmat(sigma_sqr_d + sigma_sqr_r * g .^ 2 .* abs(h_rd) .^ 2, Q, 1); % Compute the ML measurement: the weighted distance square between the received signal and the Q symbols for each of N realization
        [~, p_demod] = min(d, [], 1);
        BER(m, i_batch) = mean(B((p_demod - 1) * Q + p)) / Nbps; % BER for the m-th transmission
    end
end

BER = mean(BER, 2);