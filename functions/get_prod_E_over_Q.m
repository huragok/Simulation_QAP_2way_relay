function prod_E_over_Q_new = get_prod_E_over_Q(prod_E_over_Q, E, map)
%   prod_E_over_Q_new = get_prod_E_over_Q(prod_E_over_Q, E, map)
%   Update the product of E over Q given the update another
%   round of retransmission
% _________________________________________________________________________
%	Inputs:
%       prod_E_over_Q:      Q-by-Q matrix, the original product of E over Q
%       E:                  Q-by-Q matrix, the updating matrix in row
%                           major order.
%       map:                Q-by-Q vector, how does the Q indices are 
%                           mapped to constellation points
%	Outputs:
%       prod_E_over_Q_new:  Q-by-Q matrix, the updated product of E over Q
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 07/20/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = length(map);
prod_E_over_Q_new = zeros(Q, Q);

for p = 1 : Q
    for q = 1 : Q
        if p ~= q
            prod_E_over_Q_new(p, q) = prod_E_over_Q(p, q) * E(map(p), map(q));
        else
            prod_E_over_Q_new(p, q) = 1 / Q;
        end
    end
end