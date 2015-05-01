function BER = get_BER_upper_bound(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   BER = get_BER_upper_bound(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   Get the BER upper bounds after each retransmission based on our
%   approximation
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
%	Outputs:
%		BER:            M-by-1 vector, the BER after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/27/2015
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map); % Number of transmissions and size of the constellation
Nbps = round(log2(Q)); % Number of bit per symbol

xpcd_PBER = get_hamming_dist(Nbps) / 2 / Q / Nbps; % Initialization
BER = zeros(M, 1);
for m = 1 : M
    dist_sqr = abs(repmat(constellation(map(m, :)), 1, Q) - repmat(constellation(map(m, :)).', Q, 1)) .^ 2; % A Q-by-Q matrix containing the distance square mesurements
    dist_sqr = reshape(dist_sqr.', 1, Q^2); % Convert it to row major order
    E = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r); % Get this thing fully vectorized
    xpcd_PBER = xpcd_PBER .* E; % Update the expected number of pairwise error bit
    BER(m) = sum(sum(xpcd_PBER));
end

