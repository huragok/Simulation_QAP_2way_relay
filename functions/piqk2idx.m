function idx = piqk2idx(piqk)
%   idx = piqk2idx(piqk)
%   Convert a 4d index vector piqk into a scalar lexicographical index
% _________________________________________________________________________
%	Inputs:
%		piqk:	4-by-1 vector representing 4 indices p, i, q, k in a roll,
%               each element a value from 1 to Q
%	Outputs:
%       idx:	Scalar from 1 to Q^4, a value to index the 4D cost matrix c
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

order_inv = [2, 4, 1, 3];
piqk = piqk(:) - 1;
idx = piqk(order_inv)' * (16 .^ (3 : -1 : 0)') + 1;

