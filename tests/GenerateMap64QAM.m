clear all;
close all;
clc;

addpath('../functions');

i1 = [3, 4, 5, 6];
i2 = [2, 3, 6, 7];
i3 = [1, 2, 3, 4];

q1 = [3, 4, 5, 6];
q2 = [2, 3, 6, 7];
q3 = [5, 6, 7, 8];

loc2p_gray = [47, 46, 42, 43, 59, 58, 62, 63;
               45, 44, 40, 41, 57, 56, 60, 61;
               37, 36, 32, 33, 49, 48, 52, 53;
               39, 38, 34, 35, 51, 50, 54, 55;
               7, 6, 2, 3, 19, 18, 22, 23;
               5, 4, 0, 1, 17, 16, 20, 21;
               13, 12, 8, 9, 25, 24, 28, 29;
               15, 14, 10, 11, 27, 26, 30, 31]';
           
loc2p_new = zeros(8);

for row = 1 : 8
    for col = 1 : 8
       if (any(row == q1))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 4;
       end
       if (any(row == q2))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 2;
       end
       if (any(row == q3))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 0;
       end
       
       if (any(col == i1))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 5;
       end
       if (any(col == i2))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 3;
       end
       if (any(col == i3))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 1;
       end
    end
end

map_new = zeros(1, 64);
for idx = 0 : 63
    loc = find(loc2p_new == loc2p_gray(idx + 1));
    map_new(loc2p_gray(idx + 1) + 1) = loc2p_gray(loc) + 1;
end

%map_new = map_new + 1;

X = get_constellation(6, 'QAM', 1);
h = plot_mapping(X, map_new, 'gray', true)