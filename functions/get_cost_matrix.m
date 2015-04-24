function c = get_cost_matrix(xpcd_num_PE_bits, constellation, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   c = get_cost_matrix(xpcd_num_PE_bits, constellation, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r)
%   Compute the Q^4 4D cost matrix 
% _________________________________________________________________________
%	Inputs:
%       xpcd_num_PE_bits:   Q-BY-Q matrix, the expected number of pariwise
%                           error bits  
%       constellation:      Q-by-1 vector, the constellation of the signal 
%       beta_sr:            Scalar, the variance of the Rayleigh channel 
%                           from source to relay
%       beta_rd:            Scalar, the variance of the Rayleigh channel 
%                           from relay to destination
%       g:                  Scalar, the power normalization factor at the 
%                           relay
%       sigma_sqr_d:        Scalar, the variance of AWGN noise at the 
%                           destination      
%       sigma_sqr_r:        Scalar, the variance of AWGN noise at the relay
%	Outputs:
%		c:                  1-by-Q^4 vector, the 4D cost matrix c_piqk in
%                           the lexicalgraphical order of qpki
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________