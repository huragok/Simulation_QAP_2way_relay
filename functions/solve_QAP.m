function map = solve_QAP(a, b)
%   map = solve_QAP(f, d)
%   A QAP solver function. The QAP problem is specified in the 
%   Koopmans-Beckmann form, i.e. we would like to solve
%   \begin{align}
%       \min_x\sum_{p}\sum_{i}\sum_{q}\sum_{k}a_{pq}b_{ik}x_{pi}x_{qk} \\
%       \mbox{s.t.}\sum_{p}x_{pi} = 1,\sum_{i}x_{pi} = 1, x_{pi}\in\{0, 1\}
%   \end{align}
%   Currently just a place holder function that always returns a Gray
%   mapping.
% _________________________________________________________________________
%	Inputs:
% 		a:      1-by-Q^2 vector, the 2D "Flow" matrix a_{pq} in row major
%               order
%       b:      1-by-Q^2 vector, the 2D "Distance" matrix b_{ik} in row
%               major order
%	Outputs:
%		map:    1-by-Q vector, a permutaion of 1 : Q indicating how the Q
%               indices are mapped to constellation points
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/27/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = round(length(a) ^ (1 / 2));
map = 1 : Q;