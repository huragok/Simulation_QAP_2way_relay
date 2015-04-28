function map = solve_QAP(c)
%   map = solve_QAP(c)
%   A place holder function for the QAP solver. Always returns a Gray
%   mapping
% _________________________________________________________________________
%	Inputs:
% 		c:      1-by-Q^4 vector, the 4D cost matrix c_piqk in the
%               lexicalgraphical order of qpki
%	Outputs:
%		map:     1-by-Q vector, how does the Q indices are mapped to
%                constellation points
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 04/27/2015
% Codename: Dunkirk
% _________________________________________________________________________

Q = round(length(c) ^ (1 / 4));
map = 1 : Q;