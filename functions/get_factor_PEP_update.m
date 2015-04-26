function E = get_factor_PEP_update(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   B = get_factor_PEP_update(dist_sqr, beta, g, sigma_sqr, sigma_sqr_r)
%   Get the factor to be multiplied to PEP between indices p and q for the
%   first (M-1) retransmissions to compute the PEP for the first M 
%   retransmissions
% _________________________________________________________________________
%	Inputs:
% 		dist_sqr:       Scalar, the square norm of the 2 comstellation 
%                       points to which p and q are mapped in the M-th
%                       retransmission
%       beta_sr:        Scalart, the variance of the Rayleigh channel from
%                       source to relay
%       beta_rd:        Scalar, the variance of the Rayleigh channel from
%                       relay to destination
%       g:              Scalar, the power normalization factor at the relay
%       sigma_sqr_d:    Scalar, the variance of AWGN noise at the
%                       destination      
%       sigma_sqr_r:    Scalar, the variance of AWGN noise at the relay
%	Outputs:
%		E:              Scalar, the factor used to update PEP for the first
%                       (M-1) retransmissions into PEP for the first M
%                       retransmissions between p and q
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

E = 1

