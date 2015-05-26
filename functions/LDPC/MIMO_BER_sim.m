function BER_iter = MIMO_BER_sim(chnl_eq, max_frame, iter_max, coding_rate, nldpc, M)

% Function: simulation for MIMO system with LDPC coding and MAP detector

% example: BER_iter = MIMO_BER_sim(eye(2), 1000, 5, 3/4, 2400,2)

% chnl_eq: equivalent channel matrix G. y = G*x +n
% max_frame: number of frames in simulation 
% iter_max:maximum iteration time within the iterative receiver
% coding rate: coding rate of LDPC, {1/2,2/3,3/4,5/6}
% nldpc: bit length after channel coding, {nldpc,mod(nldpc,24)==0}
% M: Modulation type, bit per symbol, M=1: bpsk,M=2: qpsk modulation
% BER_iter: simulated bit error rate at each iteration time, a coloumn
% vector of iter_max elements
% Author: Weiliang Zeng
%modified from function LP_MIMO by Weiliang Zeng, and Mingxi Wang
%created date: 03/04/10  
% clc;clear;

 [row_num_chnl_eq,col_num_chnl_eq] =  size(chnl_eq);
 if (row_num_chnl_eq ~= 2) ||(row_num_chnl_eq ~= 2 )
     error('The dimension of chnl_eq  must be 2x2');
 end;

max_bit_error = 1000; %maximum total bit error number to stop the simulation
Ntime = 1; %1;%number of time samples at the input of linear precoder
Ns = 2; %1;%the number of transmit bitstreams
Nt = Ns;%the number of transmit antennas
Nr = 2;%the number of receiver antennas
bitpersym = M * Ntime * coding_rate;% per transmission

% file_home = pwd;
% comment = strcat('LP', int2str(M), 'point ', ' r=', num2str(coding_rate), ' L=', int2str(nldpc), ' useH=', int2str(useH), ' useG=', int2str(useG));
% f_proc = fopen(strcat(comment,'.txt'),'w');
%iter_max = 5; %maximum iteration time within the iterative receiver

%  %frequency flat fading channels
% % H = [1 0; 0 1];
% method = {'LP'};%'QAM', 'LP';

% config system
%EbNo = 10.^(SNR/10);
%EsNo = EbNo*bitpersym;% MAKE SURE!

% config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;
if coding_rate == 1
    bits_per_frame = nldpc;
else
    [H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
    bits_per_frame = length(H_cols) - length( P_matrix );
end;

% ========================================================
% config modulation
code_bit = (dec2bin([0:2^M-1])> '0') + 0;
data_mod = reshape(code_bit', 1, M*2^M);
h = modem.qammod('M',2^M, 'InputType', 'bit', 'SymbolOrder', 'gray');    % Create a modulator object
norm_pow = sqrt(mean(abs((modulate(h, data_mod'))).^2)); %factor of power normalization
sym_mod = ((modulate(h, data_mod'))/norm_pow).';%normalize the power

%generate bit vectors and mapped symbols used in MAP detector
bit_mat = (dec2bin([0:2^(Ns*M*Ntime)-1])> '0') + 0; %2^Ns*Mc*K x Ns*Mc*K
bit_mat_anti = 1- 2*bit_mat;%antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
for i = 1: Ntime*Ns
    bit_mat_tmp = reshape(bit_mat(:,(i-1)*M+1:i*M)', 1, M*2^(M*Ns*Ntime)); 
    sym_mod_mat(:,i) = modulate(h, bit_mat_tmp')/norm_pow; %2^Ns*Mc*K x Ns*K
    sym_mod_mat_real(:,i) = real(sym_mod_mat(:,i)); 
    sym_mod_mat_imag(:,i) = imag(sym_mod_mat(:,i)); 
end
sym_mod_mat = complex(sym_mod_mat_real, sym_mod_mat_imag); %convert data to complex type
% check for normalization
if (mean(abs(sym_mod).^2) - 1) > 10^4
    error('Power Normalization Error');
end;

snrpoint = 1; 
% simulation begin
% for med = 1 : length(method)
%     for snrpoint = 1:length(SNR)
%         fprintf( strcat( '\n', 'SNR', ' = %f dB\n'), SNR(snrpoint) );
%         current_time = fix(clock);
%         fprintf(  'Clock %2d:%2d:%2d\n',  current_time(4), current_time(5), current_time(6) );
%         clear error_frame;
        
%         if useG == 1
%             [val,idx_G] = min(abs(SNR(snrpoint)-[-10:5:15]));
%             G = Precoder_mat(:,:,idx_G);
%         elseif useG == 0
%             G = eye(Ntime*Ns);
%         end;
%         
%         % config Linear Precoding
%         [numM,numN] = size(G);
%         if numN ~= Ns*Ntime
%             error('The dimension of matrix G is wrong!');
%         end;
%         
        for i = 1 : max_frame
            if mod(i,5)==0
                fprintf('x');
            end;
            if mod(i,400)==0
                fprintf('\n');
            end;

            % source, coding, random interlever and modulation
            data = round( rand( Ns, bits_per_frame ) );
            if coding_rate == 1 % no coding
                for k = 1:Nt
                    codeword_inv(k) = randintrlv(data(k));
                end;
            else
                for k = 1:Ns
                    codeword(k,:) = LdpcEncode( data(k,:), H_rows, P_matrix );%ldpc encode, length == nldpc
                    codeword_inv(k,:) = randintrlv(codeword(k,:),0);% interleave   Bin_deintrlvd = randdeintrlv(Bin_intrlvd,0);
                    transmit_Mod(k,:) = ((modulate(h, codeword_inv(k,:)'))/norm_pow);
                end;
            end;

            % linear precoding
            [numM, numN] = size(transmit_Mod);
            x = reshape(transmit_Mod, Ns*Ntime, numN/Ntime);
%             s = G * x; %precoding
             [numM, numN] = size(x);
%             
            % MIMO channel
%             variance = 1/(2*EsNo(snrpoint)); %channel noise variance of I or Q path
%             H_block = kron(eye(Ntime), H); %block diagnol channel matrix with multiple time transmit vectors
            variance = 0.5;            %channel noise variance of I or Q path
            n = sqrt(variance) *complex(randn(numM, numN),randn(numM, numN)); %complex gaussian noise with unit covariance matrix
            y = chnl_eq * x + n;
            y = complex(real(y),imag(y)); %make sure y is complex, used for MAP
            
            %iterative receiver
            LextC = zeros(Nt,nldpc);
            chnl_eq = complex(real(chnl_eq),imag(chnl_eq));%convet chnl_eq to complex form of matrix
            for iter = 1 : iter_max
                %LextDemodulation_matlab = MAP_demodulate(y, LextC, chnl_eq, 2*variance, bit_mat_anti, sym_mod_mat, Nr, Ns, Ntime, M); %matlab function of MIMO MAP detector 
                LextDemodulation = MAP_demod(y, chnl_eq, bit_mat_anti, LextC, sym_mod_mat, 1.0); %C program of MIMO MAP detector, mex file. 
                for k = 1:Nt
                    LextDemo_deinv(k,:) = randdeintrlv(LextDemodulation(k,:),0); %de-interleave
                    [LLR_output_tmp errors] = MpDecode( LextDemo_deinv(k,:), H_rows, H_cols, max_iterations, decoder_type, 1, 1, data(k,:) ); %ldpc decoder
                    error_all(i,iter,k) = errors(end); %count errors in each iteration
                    LLR_output(k,:) = LLR_output_tmp(end,:);
                    LextDecoder(k,:) = LLR_output(k,:) - LextDemo_deinv(k,:); %calculate extrinic information from channel decoder
                    LextC(k,:) = randintrlv(LextDecoder(k,:), 0); %interleave
                end;
                error_frame(i,iter) = sum(error_all(i,iter,:));
            end;
            
            %count errors
            if sum(error_frame(:,end)) > max_bit_error || i == max_frame
                error_ebno(1:iter_max,snrpoint) = sum(error_frame,1)/(bits_per_frame*i*Ns);
%                 line1 = strcat('\n',int2str(SNR(snrpoint)),'\n');
%                 s = fprintf(f_proc,line1);
%                 s = fprintf(f_proc,'%f',error_ebno(end,snrpoint));
                break;
            end;
        end;
        
BER_iter = error_ebno(end,:);






