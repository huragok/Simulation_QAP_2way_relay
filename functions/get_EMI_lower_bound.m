function EMI = get_EMI_lower_bound(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   EMI = get_EMI_lower_bound(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   Get the approximated EMI lower bounds after each retransmission based 
%   on our approximation
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
%		EMI:            M-by-1 vector, the EMI after each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 07/20/2015
% Codename: Dunkirk
% _________________________________________________________________________

[M, Q] = size(map); % Number of transmissions and size of the constellation

prod_E_over_Q = 1 / Q; % Initialization
EMI = zeros(M, 1);
for m = 1 : M
    dist_sqr = abs(repmat(constellation(map(m, :)), 1, Q) - repmat(constellation(map(m, :)).', Q, 1)) .^ 2; % A Q-by-Q matrix containing the distance square mesurements
    E = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r); % Get this thing fully vectorized
    prod_E_over_Q = prod_E_over_Q .* E; % Update the expected number of pairwise error bit
    EMI(m) = log2(Q) - log2(sum(sum(prod_E_over_Q)));
end