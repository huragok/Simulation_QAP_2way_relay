function piqk = idx2piqk(idx, Q)
%   piqk = idx2piqk(idx, Q)
%   Convert a scalar lexicographical index into a 4d index vector piqk
% _________________________________________________________________________
%	Inputs: 
%       idx:	Scalar from 1 to Q^4, a value to index the 4D cost matrix c
%       Q:      Scalar, size of the constellation
%	Outputs:
%		piqk:	4-by-1 vector representing 4 indices p, i, q, k in a roll,
%               each element a value from 1 to Q
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

order = [3, 1, 4, 2];
idx_residual = idx - 1;
piqk = zeros(4, 1);
for d = 1 : 4
    piqk(order(d)) = mod(idx_residual, Q);
    idx_residual = floor(idx_residual / Q);
end
piqk = piqk + 1;