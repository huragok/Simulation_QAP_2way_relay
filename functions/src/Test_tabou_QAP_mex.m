clear all;
close all;
clc;

f = 100 * rand(3);
d = 100 * rand(3);
n_itr = 10;

[map, cost] = tabou_QAP_mex(int64(f), int64(d), int64(n_itr));