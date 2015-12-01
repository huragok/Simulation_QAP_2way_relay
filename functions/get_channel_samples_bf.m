function [h_success, g_success, h_failure, g_failure] = get_channel_samples_bf(N, constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, max_frame, iter_max, coding_rate, nldpc, seed)
%   [h_success, g_success, h_failure, g_failure] = get_channel_samples_bf(constellation, map, beta_sr, beta_rd, Pr, P1, P2, sigma_sqr_d, sigma_sqr_r, max_frame, iter_max, coding_rate, nldpc, seed)
%   Generate the samples of channels corresponding to the successful packet
%   transmissions and the failed packet transmissions, assuming block 
%   fading model.
% _________________________________________________________________________
%	Inputs:
%       N:              scalar, the number of transmissions when we
%                       decide whether a packet has been transmitted 
%                       successfully or not, must be no larger than M and
%                       larget than 0
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta_sr:        Scalar, the variance of the Rayleigh channel from
%                       source to relay
%       beta_rd:        Scalar, the variance of the Rayleigh channel from
%                       relay to destination
%       Pr:             Scalar, the average power constraint at the relay
%       P1:             Scalar, the average power constraint at the
%                       source
%       P2:             Scalar, the average power constraint at the
%                       destination
%       sigma_sqr_d:    Scalar, the variance of AWGN noise at the
%                       destination
%       sigma_sqr_r:    Scalar, the variance of AWGN noise at the relay
%       max_frame:      Scalar, number of LDPC frames in simulation.
%       iter_max:       Sclar, maximum iteration time within the iterative 
%                       receiver.
%       coding_rate:    coding rate of LDPC, {1/2,2/3,3/4,5/6}.
%       nldpc:          Scalar, bit length after channel coding, 
%                       mod(nldpc,24)=0.
%       seed:           Scalar, seed for the random number generator
%	Outputs:
%		h_success:      n_success-by-N matrix, the realization of h_1
%                       corresponding to the successful transmissions
%       g_success:      n_success-by-N matrix, the realization of g_2
%                       corresponding to the successful transmissions
%       h_failure:      n_failure-by-(N matrix, the realization of h_1
%                       corresponding to the failed transmissions
%       g_failure:      n_failure-by-N matrix, the realization of g_2
%                       corresponding to the failed transmissions
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/01/2015
% Codename: Dunkirk
% _________________________________________________________________________
% The design of this function is based on MIMO_BER_sim function by  
% Weiliang Zeng, 03/04/10.
% _________________________________________________________________________

[M, Q] = size(map);
Nbps = round(log2(Q)); % Number of bit per symbol
rng(seed);

% Config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;

[H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
bits_per_frame = length(H_cols) - length(P_matrix);

%generate bit vectors and mapped symbols used in MAP detector
bit_mat = (dec2bin(0 : Q - 1) > '0') + 0;
bit_mat_anti = 1 - 2 * bit_mat; % antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
sym_mod_mat = constellation(map(1 : N, :)).';
if (N == 1)
    sym_mod_mat = sym_mod_mat(:);
end

% Start the transmission frame by frame
h_success = zeros(0, N);
g_success = zeros(0, N);
h_failure = zeros(0, N);
g_failure = zeros(0, N);

for i = 1 : max_frame
	if mod(i, 5)==0 % Print a 'x' for every 5 frames
		fprintf('x');
	end
	if mod(i, 400)==0 % Print enter for every 400 frames
		fprintf('\n');
    end

	% source, coding, random interlever and modulation
    data = round(rand(1, bits_per_frame)); % Randomly generated information-bearing bit, 1-by-bits_per_frame
    codewordTemp = LdpcEncode(data, H_rows, P_matrix ); % LDPC encoded codeword, 1-by-nldpc
    codeword_invTemp = randintrlv(codewordTemp, 0);  % Interleaved codeword, 1-by-nldpc
    coded_index = bit2idx(codeword_invTemp, Nbps);

    transmit_Mod = zeros(N, ceil(nldpc / Nbps)); % Modulated symbol for each (re)transmission
    for m = 1 : N % 
        transmit_Mod(m, :) = constellation(map(m, coded_index)); % Modulation and normalization
    end
	
	% Lets start simulating the transmission
	numSymbol = size(transmit_Mod, 2); % Number of transmitted symbols per stream per frame
	y = zeros(N, numSymbol); % The received signal at each (re)transmission
    % Generate the channels. Assume channel to be independently fading
    % across symbols
    h_sr = sqrt(beta_sr / 2) * (randn(N, 1) + 1i * randn(N, 1)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
    h_rd = sqrt(beta_rd / 2) * (randn(N, 1) + 1i * randn(N, 1));
    g = sqrt(Pr ./ (abs(h_sr) .^ 2 * P1 + abs(h_rd) .^ 2 * P2 + sigma_sqr_r)); % The power normalization factor
    
    chnl_eq = zeros(N, numSymbol); % The equivalent channel at each (re)transmission
    for m = 1 : N 	
		% Received signal and equivalent channel generation
        y(m, :) = g(m) * h_rd(m) * (h_sr(m) * transmit_Mod(m, :) + sqrt(sigma_sqr_r / 2) * (randn(1, numSymbol) + 1i * randn(1, numSymbol))) + sqrt(sigma_sqr_d / 2) * (randn(1, numSymbol) + 1i * randn(1, numSymbol));
        
		cov = abs(g(m) * h_rd(m)) ^ 2 * sigma_sqr_r + sigma_sqr_d;
		y(m, :) = cov ^ (-1/2) * y(m, :);
		chnl_eq(m, :) = cov ^ (-1/2) * (g(m) * h_rd(m) * h_sr(m)) * ones(1, numSymbol);
        
        y(m, :) = complex(real(y(m, :)), imag(y(m, :))); % make sure y is complex, used for MAP
        chnl_eq(m, :) = complex(real(chnl_eq(m, :)), imag(chnl_eq(m, :))); % make sure chnl_eq is complex, used for MAP
    end
    
     %iterative receiver
    LextC = zeros(1, nldpc);

    error_perFrame_perMS = zeros(iter_max, 1);
    for iter = 1 : iter_max

        LextDemodulation = MAP_demod(y, chnl_eq, bit_mat_anti, LextC, sym_mod_mat, 1.0);

        LextDemo_deinv = randdeintrlv(LextDemodulation, 0); %de-interleave
        [LLR_output_tmp, errors] = MpDecode(LextDemo_deinv, H_rows, H_cols, max_iterations, decoder_type, 1, 1, data); %ldpc decoder
        error_perFrame_perMS(iter) = errors(end); %count errors in each iteration
        LLR_output = LLR_output_tmp(end, :);
        LextDecoder = LLR_output - LextDemo_deinv; %calculate extrinic information from channel decoder
        LextC = randintrlv(LextDecoder, 0); %interleave
    end
    error_all = error_perFrame_perMS(end);

    % Detection of the current frame is successful, compute the
    % effective thourhput
    if error_all == 0
        h_success = [h_success; h_sr.'];
        g_success = [g_success; h_rd.'];
    else
        h_failure = [h_failure; h_sr.'];
        g_failure = [g_failure; h_rd.'];
    end
    
end