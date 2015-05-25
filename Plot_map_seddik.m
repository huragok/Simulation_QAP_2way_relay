clear all
close all
clc

addpath('./functions/');

Nbps = 5; % 4, 5, 6
Q = 2 ^ Nbps;
M = 2;
map = get_map_seddik2(Q, M);

X = get_constellation(Nbps, 'QAM', 1);

plot_mapping(X, map(1, :), 'tmp', true);