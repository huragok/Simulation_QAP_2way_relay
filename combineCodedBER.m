clear all;
close all;
clc;

figure;
load('./waterfall_2M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'b+-', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_core, codedBER_core, 'b^-', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'bo-', 'linewidth', 2);


load('./waterfall_3M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'r+--', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_core, codedBER_core, 'r^--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'ro--', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
xlim([0, 12]), ylim([1e-5, 1])
legend({'NM1', 'CR1', 'QAP1', 'NM2', 'CR2', 'QAP2'}, 'Location', 'southwest');

figure;
load('./waterfall_4M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'b+-', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_core, codedBER_core, 'b^-', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'bo-', 'linewidth', 2);


load('./waterfall_5M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'r+--', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_core, codedBER_core, 'r^--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'ro--', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
xlim([-4, 8]), ylim([1e-5, 1])
legend({'NM3', 'CR3', 'QAP3', 'NM4', 'CR4', 'QAP4'}, 'Location', 'southwest');