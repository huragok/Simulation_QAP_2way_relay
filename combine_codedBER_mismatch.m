clear all;
close all;
clc;

cmap = [0, 0, 0; 0, 0, 1; 1, 0, 0; 0, 1, 0; 1, 1, 0];

figure;
load('./waterfall_mismatch_2M_64QAM.mat');
n_case = length(dB_inv_sigma2);
for i_case = 1 : n_case
    semilogy(dB_inv_sigma2{i_case}, codedBER{i_case}, '+-', 'Color', cmap(i_case, :), 'linewidth', 2), hold on;
end

load('./waterfall_mismatch_3M_64QAM.mat');
for i_case = 1 : n_case
    semilogy(dB_inv_sigma2{i_case}, codedBER{i_case}, '+-', 'Color', cmap(i_case, :), 'linewidth', 2), hold on;
end

grid on;
xlim([-2, 10]), ylim([1e-5, 1]);
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
legend('m=2,sol 1',...
       'm=2,sol 2',...
       'm=2,sol 3',...
       'm=2,sol 4',...
       'm=2,sol 5',...
       'm=3,sol 1',...
       'm=3,sol 2',...
       'm=3,sol 3',...
       'm=3,sol 4',...
       'm=3,sol 5')

   
figure;
load('./waterfall_mismatch_4M_64QAM.mat');
n_case = length(dB_inv_sigma2);
for i_case = 1 : n_case
    semilogy(dB_inv_sigma2{i_case}, codedBER{i_case}, '+-', 'Color', cmap(i_case, :), 'linewidth', 2), hold on;
end

load('./waterfall_mismatch_5M_64QAM.mat');
for i_case = 1 : n_case
    semilogy(dB_inv_sigma2{i_case}, codedBER{i_case}, '+-', 'Color', cmap(i_case, :), 'linewidth', 2), hold on;
end

grid on;
xlim([-4, 8]), ylim([1e-5, 1]);
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
legend('m=4,sol 1',...
       'm=4,sol 2',...
       'm=4,sol 3',...
       'm=4,sol 4',...
       'm=4,sol 5',...
       'm=5,sol 1',...
       'm=5,sol 2',...
       'm=5,sol 3',...
       'm=5,sol 4',...
       'm=5,sol 5')