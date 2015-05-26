function BER = get_codedBER(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, max_frame, iter_max, coding_rate, nldpc, seed)
%   BER = get_codedBER(constellation, map, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, max_frame, iter_max, coding_rate, nldpc)
%   Evaluate the LDPC-coded BER of a specific MoDiv mapping design over
%   multiple retransmissions.
% _________________________________________________________________________
%	Inputs:
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            M-by-Q vector, the mapping at each transmission
%       beta_sr:        Scalar, the variance of the Rayleigh channel from
%                       source to relay
%       beta_rd:        Scalar, the variance of the Rayleigh channel from
%                       relay to destination
%       g:              Scalar, the power normalization factor at the relay
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
%		BER:			M-by-1 vector, the encoded BER for each 
%                       retransmission
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 05/25/2015
% Codename: Dunkirk
% _________________________________________________________________________
% The design of this function is based on MIMO_BER_sim function by  
% Weiliang Zeng, 03/04/10.
% _________________________________________________________________________

[M, Q] = size(map);
Nbps = round(log2(Q)); % Number of bit per symbol
symbols_mapped = constellation(map); % The mapped symbols at all transmissions

max_bit_error = 1000; % Count up to this number of bit error we stop the simulation since BER can be measured accurately enough at this point

% Config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;

[H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
bits_per_frame = length(H_cols) - length(P_matrix);

%generate bit vectors and mapped symbols used in MAP detector
bit_mat = (dec2bin(0 : Q - 1) > '0') + 0;
sym_mod_mat = constellation(map);
for i = 1 : d
    bit_mat_tmp = reshape(bit_mat(:,(i - 1) * Nbps + 1 : i * Nbps)', 1, Nbps * 2 ^ (Nbps * d)); 
    sym_mod_mat(:, i) = modulate(h, bit_mat_tmp') / norm_pow; %2^Ns*Mc*K x Ns*K
end

% Generate the channels. Assume a block fading channel which is  stationary
% within each frame and independently fading across frames
h_sr = sqrt(beta_sr / 2) * (randn(max_frame, M) + 1i * randn(max_frame, M)); % Generate the Rayleigh channel, We expect N to be large so inorder to cut memory usage we generate the random channel/noise once for each transmission
h_rd = sqrt(beta_rd / 2) * (randn(max_frame, M) + 1i * randn(max_frame, M));
        
% Start the transmission frame by frame
error_all = zeros(M, max_frame);
for i = 1 : max_frame
	if mod(i, 5)==0 % Print a 'x' for every 5 frames
		fprintf('x');
	end
	if mod(i, 400)==0 % Print enter for every 400 frames
		fprintf('\n');
	end
	
	% source, coding, random interlever and modulation

    data = round(rand(1, bits_per_frame)); % Randomly generated information-bearing bit, d-by-bits_per_frame
    codewordTemp = LdpcEncode(data, H_rows, P_matrix ); % LDPC encoded codeword, length = nldpc
    codeword_invTemp = randintrlv(codewordTemp, 0);  % Interleaved codeword, length = nldpc
    coded_index = bit2idx(codeword_invTemp);

    transmit_Mod = zeros(M, ceil(nldpc / Nbps)); % Modulated symbol for each (re)transmission
    for m = 1 : M % 
        transmit_Mod(m, :) = constellation(map(m, coded_index)); % Modulation and normalization
    end
	
	% Lets start simulating the transmission
	numSymbol = size(transmit_Mod, 2); % Number of transmitted symbols per stream per frame
	y = zeros(M, numSymbol); % The received signal at each (re)transmission
    chnl_eq = zeros(M, 1); % The equivalent channel at each (re)transmission
    for m = 1 : M	
		% Received signal and equivalent channel generation
        y(m, :) = g * h_rd(i, m) * (h_sr(i, m) * transmit_Mod(m, :) + sqrt(sigma_sqr_r / 2) * (randn(1, numSymbol) + 1i * randn(1, numSymbol))) + sqrt(sigma_sqr_d / 2) * (randn(1, numSymbol) + 1i * randn(1, numSymbol));
        
		cov = abs(g * h_rd(i, m)) ^ 2 * sigma_sqr_r + sigma_sqr_d;
		y(m, :) = cov ^ (-1/2) * y(m, :);
		chnl_eq(m) = cov ^ (-1/2) * g * h_rd(i, m) * h_sr(i, m);
    end
    y = complex(real(y), imag(y)); % make sure y is complex, used for MAP
	chnl_eq = complex(real(chnl_eq), imag(chnl_eq)); % make sure chnl_eq is complex, used for MAP
    
    for m = 1 : M	
		%iterative receiver
        LextC = zeros(1, nldpc);
		
		error_perFrame_perMS = zeros(iter_max, 1);
		for iter = 1 : iter_max

            LextDemodulation = MAP_demod(y(1 : m, :), chnl_eq(1 : m), bit_mat_anti, LextC, sym_mod_mat, 1.0);

            LextDemo_deinv = randdeintrlv(LextDemodulation, 0); %de-interleave
            [LLR_output_tmp, errors] = MpDecode(LextDemo_deinv, H_rows, H_cols, max_iterations, decoder_type, 1, 1, data); %ldpc decoder
            error_perFrame_perMS(iter) = errors(end); %count errors in each iteration
            LLR_output = LLR_output_tmp(end, :);
            LextDecoder = LLR_output - LextDemo_deinv; %calculate extrinic information from channel decoder
            LextC = randintrlv(LextDecoder, 0); %interleave
		end

		error_all(m, i) = error_perFrame_perMS(end);  
    end
    %count errors
    if min(sum(error_all, 2)) > max_bit_error
        break;
    end
end

BER = sum(error_all, 2) / (bits_per_frame * i);
