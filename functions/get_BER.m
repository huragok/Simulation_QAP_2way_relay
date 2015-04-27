function BER = get_BER(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, N)
%   BER = get_BER(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, N)
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
%       g:              Scalar, the power normalization factor at the relay
%       sigma_sqr_d:    Scalar, the variance of AWGN noise at the
%                       destination
%       sigma_sqr_r:    Scalar, the variance of AWGN noise at the relay
%       N:              Scalar, number of Monte-Carlo run
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
p = randi(Q, 1, N); % Generate the random transmitted index
d = zeros(Q, N); % The distance measurement used for the ML demodulator
B = get_hamming_dist(Nbps); % The hamming distance matrix

BER = zeros(M, 1);
for m = 1 : M
    h_sr = sqrt(beta_sr / 2) * (randn(1, N) + 1i * randn(1, N)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
    h_rd = sqrt(beta_rd / 2) * (randn(1, N) + 1i * randn(1, N));
    symbol = symbols_mapped(m, p); % Generate the random symbol
    
    y = g * h_rd .* (h_sr .* symbol + sqrt(sigma_sqr_r / 2) * (randn(1, N) + 1i * randn(1, N))) + sqrt(sigma_sqr_d / 2) * (randn(1, N) + 1i * randn(1, N));
    
    d = d + abs(repmat(y, Q, 1) - g * symbols_mapped(m, :).' * (h_rd .* h_sr)) .^ 2 ./ repmat(sigma_sqr_d + g ^ 2 * sigma_sqr_r * abs(h_rd) .^ 2, Q, 1); % Compute the ML measurement: the weighted distance square between the received signal and the Q symbols for each of N realization
    [~, p_demod] = min(d, [], 1);
    idx_B = (p_demod - 1) * Q + p; % The 1-d vector index for B corresponding to the 2d index (p, p_demod)
    BER(m) = mean(B(idx_B)) / Nbps; % BER for the m-th transmission
end