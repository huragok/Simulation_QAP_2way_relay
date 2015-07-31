clear all;
close all;
clc;

addpath('../functions');

i1 = [1, 4];
i2 = [1, 2];

q1 = [1, 4];
q2 = [3, 4];

loc2p_gray = [11, 10, 14, 15, 9, 8, 12, 13, 1, 0, 4, 5, 3, 2, 6, 7]; % In column major order
           
loc2p_new = zeros(4);

for row = 1 : 4
    for col = 1 : 4
       if (any(row == q1))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 2;
       end
       if (any(row == q2))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 0;
       end

       
       if (any(col == i1))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 3;
       end
       if (any(col == i2))
           loc2p_new(row, col) = loc2p_new(row, col) + 2 ^ 1;
       end
    end
end

map_new = zeros(1, 16);
for idx = 0 : 15
    loc = find(loc2p_new == loc2p_gray(idx + 1));
    map_new(loc2p_gray(idx + 1) + 1) = loc2p_gray(loc) + 1;
end

%map_new = map_new + 1;

X = get_constellation(4, 'QAM', 1);
h = plot_mapping(X, map_new, 'gray', true)