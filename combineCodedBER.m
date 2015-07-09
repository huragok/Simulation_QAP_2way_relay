clear all;
close all;
clc;

figure;
load('waterfall_3M_16QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'b+-', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_seddik, codedBER_seddik, 'b^-', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'bo-', 'linewidth', 2);

load('waterfall_4M_16QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'r+--', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_seddik, codedBER_seddik, 'r^--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'ro--', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
xlim([-4, 8]), ylim([1e-5, 1])
legend({'NM3', 'SE3', 'QAP3', 'NM3', 'SE4', 'QAP4'}, 'Location', 'southwest');