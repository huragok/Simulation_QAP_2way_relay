function [sym_mod, sym_mod_mat] = modConfig(Ch)

% config modulation
code_bit = (dec2bin([0:2^Ch.M-1])> '0') + 0;
data_mod = reshape(code_bit', 1, Ch.M*2^Ch.M);
h = modem.qammod('M',2^Ch.M, 'InputType', 'bit', 'SymbolOrder', 'gray');    % Create a modulator object
norm_pow = sqrt(mean(abs((modulate(h, data_mod'))).^2)); %factor of power normalization
sym_mod = ((modulate(h, data_mod'))/norm_pow).';%normalize the power

%generate bit vectors and mapped symbols used in MIMO detector
bit_mat = (dec2bin([0:2^(Ch.Ns*Ch.M)-1])> '0') + 0; %(2^Ns*Ch.M*Ch.Ntone) x (Ns*Ch.M*Ch.Ntone),
bit_mat_anti = 1- 2*bit_mat;%antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
for i = 1: Ch.Ns
	bit_mat_tmp = reshape(bit_mat(:,(i-1)*Ch.M+1:i*Ch.M)', 1, Ch.M*2^(Ch.M*Ch.Ns));
	sym_mod_mat(:,i) = modulate(h, bit_mat_tmp')/norm_pow; %2^Ch.Ns*Mc*K x Ch.Ns*K
end

if (mean(abs(sym_mod).^2) - 1) > 10^4
	error('Power Normalization Error');
end;