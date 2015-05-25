clear all;
close all;
clc;

along_gray_MSB = [0, 1, 3, 2, 6, 7, 5, 4];
along_seddik_MSB = [5, 0, 6, 3, 2, 7, 4, 1];

along_gray_LSB = [0, 1, 3, 2];
along_seddik_LSB = [1, 2, 0, 3];

const_gray = zeros(4, 8);
const_seddik = zeros(4, 8);
for i = 1 : 4
    for j = 1 : 8
        const_gray(i, j) = along_gray_MSB(j) * 4 + along_gray_LSB(i);
        const_seddik(i, j) = along_seddik_MSB(j) * 4+ along_seddik_LSB(i);
    end
end

map_seddik = zeros(1, 32);
for idx = 0 : 31
    loc = find(const_seddik == idx);
    map_seddik(idx + 1) = const_gray(loc);
end

map_seddik = map_seddik + 1;