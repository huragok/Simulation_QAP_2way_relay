clear all;
close all;
clc;

addpath('functions');

%% 1. Simulation settings
% 64QAM
load('./data/Test_201573023658.mat');
dB_inv_sigma2 = 13.5; % 13.5, 5.8, 2, -0.2
N = 1; % 1, 2, 3, 4

Nbps = test_cases(1).param_origin.Nbps;
type_mod = test_cases(1).param_origin.type_mod;
p_Pr = test_cases(1).param_origin.p_Pr; % Power at the relay
d = [test_cases(1).param_origin.d1, test_cases(1).param_origin.d2]; % Distance between S and R, R and D
nu = test_cases(1).param_origin.nu; % Pathloss factor
M = test_cases(1).param_origin.M; % Total number of transmissions

max_frame = 20;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

%% 2. Initialization
Q = 2 ^ Nbps;
Pr = 4 * p_Pr;
Pt = 2 * (1 - p_Pr);
constellation = sqrt(Pt) * get_constellation(Nbps, type_mod, 1);

beta_sr = d(1) ^ -nu;
beta_rd = d(2) ^ -nu;

map  = [1 : Q; test_cases(1).map];

sigma2 = 10 .^ (-dB_inv_sigma2 / 10); % The noise covariance at all nodes

%% 3. Generate the channel samples corresponding to the successful transmission and the failed transmissions
[h_success, g_success, h_failure, g_failure] = get_channel_samples(N, constellation, map, beta_sr, beta_rd, Pr, Pt, Pt, sigma2, sigma2, max_frame, iter_max, coding_rate, nldpc, seed);

%% 4. Visualization, plot the mean and covaraince matrix
points_success = [abs(h_success), abs(g_success)];
points_failure = [abs(h_failure), abs(g_failure)];

% The mean plot
mean_success = mean(points_success);
mean_failure = mean(points_failure);
figure;
plot(1 : 2 * N, mean_success, 'bo-', 'linewidth', 2), hold on;
plot(1 : 2 * N, mean_failure, 'r+--', 'linewidth', 2);
grid on;
set(gca, 'Fontsize', 18);
set(gca, 'XTick', 1 : 2 * N);
xlabel('[|h(1 : N)|, |g(1, N)|]'), ylabel('mean');
ylim([0, ceil(max(max(mean_success), max(mean_failure)))]);
legend('Success', 'Failure');

% The covariance plot
colormap('hot');

cov_success = cov(points_success);
cov_failure = cov(points_failure);

figure;
imagesc(cov_success);
axis equal;
xlabel('[|h(1 : N)|, |g(1, N)|]'), ylabel('[|h(1 : N)|, |g(1, N)|]');
xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
set(gca, 'XTick', 1 : 2 * N);
set(gca, 'YTick', 1 : 2 * N);
set(gca, 'Fontsize', 18);
colorbar;
caxis([0, max(max(max(cov_success)), max(max(cov_failure)))]);

figure;
imagesc(cov_failure);
axis equal;
xlabel('[|h(1 : N)|, |g(1, N)|]'), ylabel('[|h(1 : N)|, |g(1, N)|]');
xlim([0.5, 2 * N + 0.5]), ylim([0.5, 2 * N + 0.5]);
set(gca, 'XTick', 1 : 2 * N);
set(gca, 'YTick', 1 : 2 * N);
set(gca, 'Fontsize', 18);
colorbar;
caxis([0, max(max(max(cov_failure)), max(max(cov_failure)))]);

%% 5. Output the result to a xml file used for normality test, each row corresponds to 1 observation in the order of real(h(1)), imag(h(1)), ...real(h(N)), imag(h(N)), real(g(1)), imag(g(1)), ...real(g(N)), imag(g(N)) separated by comma
doc = com.mathworks.xml.XMLUtils.createDocument('ChannelSamples');

channel_samples = doc.getDocumentElement;
channel_samples.setAttribute('nSuccess', num2str(size(h_success, 1)));
channel_samples.setAttribute('nFailure', num2str(size(h_failure, 1)));
channel_samples.setAttribute('beta', num2str(beta_sr));
channel_samples.setAttribute('N', num2str(N));

% The success part
success = doc.createElement('Success');
for i_success = 1 : size(h_success, 1)
    entry = doc.createElement('Entry');
    entry.appendChild(doc.createTextNode(num2str([real(h_success(i_success, :)), imag(h_success(i_success, :)), real(g_success(i_success, :)), imag(g_success(i_success, :))], '%-f ')));
    success.appendChild(entry);
end
channel_samples.appendChild(success);

% The failure part
failure = doc.createElement('failure');
for i_failure = 1 : size(h_failure, 1)
    entry = doc.createElement('Entry');
    entry.appendChild(doc.createTextNode(num2str([real(h_failure(i_failure, :)), imag(h_failure(i_failure, :)), real(g_failure(i_failure, :)), imag(g_failure(i_failure, :))], '%-f ')));
    failure.appendChild(entry);
end
channel_samples.appendChild(failure);

xmlwrite('samples.xml', doc);
%% 6. In the low dimensional case, a 2-D cluster visualization
if N == 1
    figure;
    scatter(points_success(:, 1), points_success(:, 2), 'r^'), hold on;
    scatter(points_failure(:, 1), points_failure(:, 2), 'bo');
    grid on;
    axis equal;
    set(gca, 'Fontsize', 18);
    xlabel('|h|'), ylabel('|g|');
    legend({'Success', 'Failure'}, 'Location', 'northeast');
end