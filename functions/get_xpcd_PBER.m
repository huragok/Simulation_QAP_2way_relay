function xpcd_PBER = get_xpcd_PBER(c, map)
%   xpcd_PBER = get_xpcd_PBER(c, map)
%   Compute the expected pairwise bit error rate given the cost matrix and
%   mapping scheme
% _________________________________________________________________________
%	Inputs:
%       c:                  1-by-Q^4 vector, the 4D cost matrix c_piqk in
%                           the lexicalgraphical order of qpki
%       map:                1-by-Q vector, how does the Q indices are 
%                           mapped to constellation points
%	Outputs:
%       xpcd_PBER:          Q-BY-Q matrix, the expected pairwise bit error
%                           rate
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/24/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = length(map);
xpcd_PBER = zeros(Q, Q);

for p = 1 : Q
    for q = 1 : Q
        if p ~= q
            piqk = [p, map(p), q, map(q)];
            xpcd_PBER(p, q) = c(piqk2idx(piqk));
        end
    end
end