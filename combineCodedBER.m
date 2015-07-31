clear all;
close all;
clc;

load('waterfall_EMI_3M_64QAM.mat');
dB_inv_sigma2_QAP_new = dB_inv_sigma2_QAP;
codedBER_QAP_new = codedBER_QAP;

figure;
load('./data/waterfall_3M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'b+-', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_seddik, codedBER_seddik, 'b^-', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'bo-', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP_new, codedBER_QAP_new, 'bs-', 'linewidth', 2);

load('waterfall_EMI_4M_64QAM.mat');
dB_inv_sigma2_QAP_new = dB_inv_sigma2_QAP;
codedBER_QAP_new = codedBER_QAP;

load('./data/waterfall_4M_64QAM.mat');
semilogy(dB_inv_sigma2_noncore, codedBER_noncore, 'r+--', 'linewidth', 2), hold on;
semilogy(dB_inv_sigma2_seddik, codedBER_seddik, 'r^--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP, codedBER_QAP, 'ro--', 'linewidth', 2);
semilogy(dB_inv_sigma2_QAP_new, codedBER_QAP_new, 'rs--', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('1/\sigma^2(dB)'), ylabel('Coded BER');
xlim([-4, 8]), ylim([1e-5, 1])
legend({'NM2', 'GS2', 'QAP2-BER', 'QAP2-EMI', 'NM3', 'GS3', 'QAP3-BER', 'QAP3-EMI'}, 'Location', 'southwest');