function bound = get_GLB(f, d)
%   bound = get_GLB(f, d)
%   Evaluate the Gilmore-Lawler bound for KB-form QAP
% _________________________________________________________________________
%	Inputs:
%       f:      Q-by-Q non-negative matrix, the 'flow' matrix
%       d:      Q-by-Q non-negative matrix, the 'distance' matrix
%
%	Outputs:
%       bound:  Scalar, the Gilmore-Lawler bound
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 09/15/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = size(f, 1);

f_hat = zeros(Q, Q - 1); % The q-th row is f(q, :) excluding f(q, q) sorted in ascending order
for q = 1 : Q
    f_hat(q, :) = sort([f(q, 1 : q - 1), f(q, q + 1 : Q)]);
end

d_hat = zeros(Q, Q - 1); % The q-th row is d(q, :) excluding d(q, q) sorted in descending order
for q = 1 : Q
    d_hat(q, :) = sort([d(q, 1 : q - 1), d(q, q + 1 : Q)], 'descend');
end

l = f_hat * d_hat'; % The matrix used for the linear assignment problem that will result in GLB

[~, bound] = munkres(l); % Hungarian algorithm