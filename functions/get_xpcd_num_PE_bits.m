function xpcd_num_PE_bits = get_xpcd_num_PE_bits(c, map)
%   xpcd_num_PE_bits = get_xpcd_num_PE_bits(c, map)
%   Compute the expected number of pairwise error bits given the cost
%   matrix and mapping scheme
% _________________________________________________________________________
%	Inputs:
%       c:                  1-by-Q^4 vector, the 4D cost matrix c_piqk in
%                           the lexicalgraphical order of qpki
%       map:                1-by-Q vector, how does the Q indices are 
%                           mapped to constellation points
%	Outputs:
%       xpcd_num_PE_bits:   Q-BY-Q matrix, the expected number of pariwise
%                           error bits  
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________