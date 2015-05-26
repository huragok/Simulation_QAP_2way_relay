function BER = MCMIMO_BER_sim(H, P, sigma2, optBD, max_frame, iter_max, coding_rate, nldpc, Nbps)

% Filename: TestBER.m
% File description:
%   Evaluate the BER of a specific precoder design algorithm over
%   a range of noise variance to a txt file in our multi-cell MIMO
%   settings.
% _________________________________________________________________________
% Output arguments:
%   BER:        	Scalar, the bit error rate
% Input arguments:
%   H:          	K-by-1 cell of Nr-by-NtTotal matrix, H{k} is the 
%					channel matrix to the k-th MS.
%   P:          	K-by-1 cell of NtTotal-by-dt matrix, P{k} is the 
%					precoder for the k-th MS.
%   sigma2:     	Scalar, the noise variance.
%   optBD:      	Bool scalar, if optBD we verify the ZF of the precoders
%					If ~optBD the interference plus noise is whitened to 
%					get the equivalent channel.
% 	max_frame: 		Scalar, number of LDPC frames in simulation.
% 	iter_max:		Sclar, maximum iteration time within the iterative 
%					receiver.
% 	coding_rate: 	coding rate of LDPC, {1/2,2/3,3/4,5/6}.
% 	nldpc: 			bit length after channel coding, mod(nldpc,24)=0.
%   Nbps:       	scalar, number of bit per symbol, determines which QAM
%               	constellation to use.
% _________________________________________________________________________
% Current version: v1.0
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 02/21/2014
% Codename: Ritornello
% _________________________________________________________________________
% The design of this function is based on MIMO_BER_sim function by  
% Weiliang Zeng, 03/04/10.
% _________________________________________________________________________

K = length(H);
[Nr, NtTotal] = size(H{1});
if (Nr ~= size(P{1}, 2) || NtTotal ~= size(P{1}, 1))
    error('The channel and precoder size dont match!');
end
if (Nr ~= 2)
    error('Sorry this program only support Nr = 2!');
end

max_bit_error = 1000; %maximum total bit error number to stop the simulation
d = Nr; % the number of transmit bitstreams

% The interference covariance from the l-th transmission to the k-th MS
covLeak = cell(K, K);
pwrLeak = zeros(K, K);
for k = 1 : K
	for l = 1 : K
		if l ~= k
			HPTemp = H{k} * P{l};
			covLeak{k, l} = (HPTemp * HPTemp');
			pwrLeak(k, l) = trace(covLeak{k, l});
		end
	end
end
if optBD
	disp(['maximum leakage power is: ', num2str(max(max(pwrLeak)))]);
end


% Config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;

[H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
bits_per_frame = length(H_cols) - length(P_matrix);

% config modulation
code_bit = (dec2bin([0 : 2 ^ Nbps - 1]) > '0') + 0;
data_mod = reshape(code_bit', 1, Nbps * 2 ^ Nbps);
h = modem.qammod('M', 2 ^ Nbps, 'InputType', 'bit', 'SymbolOrder', 'gray'); % Create a modulator object
norm_pow = sqrt(mean(abs((modulate(h, data_mod'))) .^ 2)); % Factor of power normalization
sym_mod = ((modulate(h, data_mod')) / norm_pow).'; % Normalize the power so that the we have E[ss'] = I

%generate bit vectors and mapped symbols used in MAP detecto
bit_mat = (dec2bin([0 : 2 ^ (d * Nbps) - 1]) > '0') + 0;
bit_mat_anti = 1 - 2 * bit_mat; % antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
for i = 1 : d
    bit_mat_tmp = reshape(bit_mat(:,(i - 1) * Nbps + 1 : i * Nbps)', 1, Nbps * 2 ^ (Nbps * d)); 
    sym_mod_mat(:, i) = modulate(h, bit_mat_tmp') / norm_pow; %2^Ns*Mc*K x Ns*K
end
sym_mod_mat = complex(real(sym_mod_mat), imag(sym_mod_mat)); % Make sure that the modulated symbol is comlex
% check for normalization
if (mean(abs(sym_mod).^2) - 1) > 10^4
    error('Power Normalization Error');
end

% Start the transmission frame by frame
error_all = zeros(max_frame, K * d);
for i = 1 : max_frame
	if mod(i, 5)==0 % Print a 'x' for every 5 frames
		fprintf('x');
	end
	if mod(i, 400)==0 % Print enter for every 400 frames
		fprintf('\n');
	end
	
	% source, coding, random interlever and modulation
	data = cell(K, 1); % The information bearing bits for each MS
	transmit_Sum = zeros(NtTotal, ceil(nldpc / Nbps)); % The LDPC-coded, interleaved, modulated, precoded and summed transmitted signal
	for k = 1 : K
		data{k} = round(rand(d, bits_per_frame)); % Randomly generated information-bearing bit, d-by-bits_per_frame
		transmit_Mod = zeros(d, ceil(nldpc / Nbps));
		for indStream = 1 : d
			codewordTemp = LdpcEncode(data{k}(indStream, :), H_rows, P_matrix ); % LDPC encoded codeword, length = nldpc
			codeword_invTemp = randintrlv(codewordTemp, 0);  % Interleaved codeword, length = nldpc
			transmit_Mod(indStream, :) = ((modulate(h, codeword_invTemp')) / norm_pow); % Modulation and normalization
		end
		transmit_Sum = transmit_Sum + P{k} * transmit_Mod;
	end
	
	% Lets start simulating the transmission
	numSymbol = size(transmit_Sum, 2); % Number of transmitted symbols per stream per frame
	for k = 1 : K	
		% Received signal and equivalent channel generation
		noise = sqrt(sigma2 / 2) * complex(randn(d, numSymbol), randn(d, numSymbol));
		y = H{k} * transmit_Sum + noise; % The received signal at the k-th MS
		if optBD
			chnl_eq = H{k} * P{k}; % Equivalent channel is just the multiplication of the 
		else % Here we assume that the receiver whiten the interference + noise then demulate
			cov = zeros(Nr, Nr);
			for l = 1 : K
				if l ~= k
					cov = cov + covLeak{k, l};
				end
			end
			cov = cov + sigma2 * eye(Nr);
			y = cov ^ (-1/2) * y;
			chnl_eq = cov ^ (-1/2) * H{k} * P{k};
		end
		y = complex(real(y), imag(y)); %make sure y is complex, used for MAP
		chnl_eq = complex(real(chnl_eq), imag(chnl_eq)); %make sure chnl_eq is complex, used for MAP
		
		%iterative receiver
        LextC = zeros(d, nldpc);
		
		error_perFrame_perMS = zeros(iter_max, d);
		for iter = 1 : iter_max
            if optBD
                LextDemodulation = MAP_demod(y, chnl_eq, bit_mat_anti, LextC, sym_mod_mat, sigma2); %C program of MIMO MAP detector, mex file. 
            else
                LextDemodulation = MAP_demod(y, chnl_eq, bit_mat_anti, LextC, sym_mod_mat, 1.0);
            end
            for indStream = 1 : d
				LextDemo_deinv(indStream, :) = randdeintrlv(LextDemodulation(indStream, :), 0); %de-interleave
				[LLR_output_tmp errors] = MpDecode(LextDemo_deinv(indStream, :), H_rows, H_cols, max_iterations, decoder_type, 1, 1, data{k}(indStream, :)); %ldpc decoder
				error_perFrame_perMS(iter, indStream) = errors(end); %count errors in each iteration
				LLR_output(indStream, :) = LLR_output_tmp(end, :);
				LextDecoder(indStream, :) = LLR_output(indStream, :) - LextDemo_deinv(indStream, :); %calculate extrinic information from channel decoder
				LextC(indStream, :) = randintrlv(LextDecoder(indStream, :), 0); %interleave
			end
		end
		for indStream = 1 : d
			error_all(i, (k - 1) * d + indStream) = error_perFrame_perMS(end, indStream);
		end
	end
	
	%count errors
	if sum(sum(error_all)) > max_bit_error
		break;
	end
end

% sum(error_all, 1)
% figure;
% stem3(pwrLeak);

BER = sum(sum(error_all)) / (bits_per_frame * i * d * K);
