clear all;
close all;
clc;

along_gray = [0, 1, 3, 2, 6, 7, 5, 4];
along_seddik = [2, 5, 3, 0, 4, 7, 1, 6];

const_gray = zeros(8);
const_seddik = zeros(8);
for i = 1 : 8
    for j = 1 : 8
        const_gray(i, j) = along_gray(j) * 8 + along_gray(i);
        const_seddik(i, j) = along_seddik(j) * 8 + along_seddik(i);
    end
end

map_seddik = zeros(1, 64);
for idx = 0 : 63
    loc = find(const_seddik == idx);
    map_seddik(idx + 1) = const_gray(loc);
end

map_seddik = map_seddik + 1;