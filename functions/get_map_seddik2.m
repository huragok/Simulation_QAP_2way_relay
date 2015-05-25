function map = get_map_seddik2(Q, M)
%   map = get_map_noncore(Q, M)
%   Get the Seddik's remapping for a Q-QAM constellation: i.e. use Gray
%   mapping for the first transmission and use Seddik's remapping for the
%   M-1 retransmissions
% _________________________________________________________________________
%	Inputs:
% 		Q:      scalar, size of the constellation Q must be 16, 32 or 64
%       M:      scalar, total number of transmissions
%	Outputs:
%		map:    M-by-Q matrix, the mapping at each transmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/14/2015
% Codename: Dunkirk
% _________________________________________________________________________

seddik_16 = [6, 7, 8, 5, 10, 11, 12, 9, 14, 15, 16, 13, 2, 3, 4, 1];
%seddik_32 = [12, 9, 10, 11, 24, 21, 22, 23, 4, 1, 2, 3, 16, 13, 14, 15, 28, 25, 26, 27, 8, 5, 6, 7, 20, 17, 18, 19, 32, 29, 30, 31];
%seddik_32 = [8,5,6,7,20,17,18,19,28,25,26,27,12,9,10,11,24,21,22,23,4,1,2,3,16,13,14,15,32,29,30,31];
seddik_32 = [8,5,6,7,20,17,18,19,32,29,30,31,12,9,10,11,24,21,22,23,4,1,2,3,16,13,14,15,28,25,26,27];
seddik_64 = [19, 22, 17, 20, 23, 18, 21, 24, 43, 46, 41, 44, 47, 42, 45, 48, 3, 6, 1, 4, 7, 2, 5, 8, 27, 30, 25, 28, 31, 26, 29, 32, 51, 54, 49, 52, 55, 50, 53, 56, 11, 14, 9, 12, 15, 10, 13, 16, 35, 38, 33, 36, 39, 34, 37, 40, 59, 62, 57, 60, 63, 58, 61, 64];

if (Q == 16)
    seddik = seddik_16;
elseif (Q == 32)
    seddik = seddik_32;
elseif (Q == 64)
    seddik = seddik_64;
else
    error('Constellation size must be 16, 32 or 64.');
end

map = [1 : Q; repmat(seddik, M - 1, 1)];