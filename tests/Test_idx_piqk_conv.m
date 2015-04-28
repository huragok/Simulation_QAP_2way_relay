% Function to test the conversion between piqk and idx

clear all;
close all;
clc;

addpath('../functions/');

Q = 16;
piqk = [1,13,2,16];
idx = piqk2idx(piqk)
piqk_new = idx2piqk(idx, Q)