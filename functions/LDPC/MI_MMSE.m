function [I_finite, MMSE] = MI_MMSE(Ch, eqH)
    
% this function is used for calculate the mutual information for finite
% input distribution via integral.
    
% configuration for integration
if Ch.MI == 1
    if Ch.M == 1
        Nint = 15;% 21; %35 seconds  %25; %31; %35 number of integration points for each component
    elseif Ch.M == 2
        Nint = 12; %13; %65 seconds  %15; %102 seconds
    elseif Ch.M == 3
        Nint=9; %402 seconds  %9; %294 seconds  %11; %543 seconds
    elseif Ch.M == 4
        Nint=9;
    elseif Ch.M == 6
        Nint=8;
    else
        Nint=7;
    end
    lim = 6;           		%integration limit
elseif Ch.MMSE == 1
    if Ch.M == 1
        Nint = 21;% 21; %35 seconds  %25; %31; %35 number of integration points for each component
    elseif Ch.M == 2
        Nint = 17; %13; %65 seconds  %15; %102 seconds
    elseif Ch.M == 3
        Nint=9; %402 seconds  %9; %294 seconds  %11; %543 seconds
    elseif Ch.M == 4
        Nint=9;
    elseif Ch.M == 6
        Nint=8;
    else
        Nint=7;
    end
    lim = 8;           		%integration limit
end;

delta= 2 * lim / (Nint-1);      %integration step size
smp = [-lim:delta:lim];		%integration sampling points
n_smp = length(smp);
[noiseR, noiseI] = meshgrid(smp, smp);
noise = reshape(noiseR + 1j * noiseI, 1, n_smp^2);
[noise_a, noise_b] = meshgrid(noise, noise);
noise_vec = [reshape(noise_a, 1, n_smp^4); reshape(noise_b, 1, n_smp^4)];
x_mod = Ch.sym_mod_mat.';

if Ch.MI == 1
	for m = 1:2^(2*Ch.M) % for the outer loop
		f_inner = 0;
		for k = 1:2^(2*Ch.M)
			x_diff(:,k) = x_mod(:,m)-x_mod(:,k);
			f_inner = f_inner + exp(- sum(abs(eqH*x_diff(:,k)*ones(1,n_smp^4)+noise_vec).^2, 1) ...
				+ sum(abs(noise_vec).^2,1) );
		end;
		F_avg_n(m) = sum( log2(f_inner) .* exp(-sum(abs(noise_vec).^2,1))/pi^2 )*delta^4;
	end;
	I_finite = 2*log2(2*Ch.M) - 1/(2^(2*Ch.M)) * sum(F_avg_n);
else
	I_finite = 0;
end;

if Ch.MMSE == 1
	% version 1: using matrix calculation the fastest one.
	clear f_inner
	f_avg = 0;
    f_avg_l = 0;
	for l = 1:2^(2*Ch.M)
		f_inner1 = 0;
		f_inner2 = 0;
		f_inner3 = 0;
		for k = 1:2^(2*Ch.M); % M = 2^Ch.M => M^2 = (2^Ch.M)^2 = 2^(Ch.M*2)
			f_inner1 = f_inner1 + x_mod(:,k)*exp(- ...
				(sum(abs(eqH*(x_mod(:,l)-x_mod(:,k)) *ones(1,n_smp^4) + noise_vec).^2,1) ));
		end;
		f_inner2 = f_inner1';
		
		for m = 1:2^(2*Ch.M)
			f_inner3 = f_inner3 + exp(- ...
				(sum(abs(eqH*(x_mod(:,l)-x_mod(:,m)) *ones(1,n_smp^4) + noise_vec).^2,1) ));
		end;
		
        f_inner1 = f_inner1.*(ones(2,1) * (exp(-sum(abs(noise_vec).^2,1))/pi^2 ./ (f_inner3).^2));
		f_avg_l = f_avg_l + f_inner1 * f_inner2 * delta^4;
	end;
	MMSE = eye(2) - 1/(2^(2*Ch.M)) * f_avg_l;
	
	% version 2: int first sum second (much faster than verison 3) : 2.386844s
elseif Ch.MMSE == 2
	clear f_inner
	f_avg = 0;
	for l = 1:2^(2*Ch.M)
		f_inner1 = 0;
		f_inner2 = 0;
		f_inner3 = 0;
		for k = 1:2^(2*Ch.M)
			f_inner1 = f_inner1 + x_mod(:,k)*exp(- ...
				(sum(abs(eqH*(x_mod(:,l)-x_mod(:,k)) *ones(1,n_smp^4) + noise_vec).^2,1) ));
		end;
		f_inner2 = f_inner1';
		
		for m = 1:2^(2*Ch.M)
			f_inner3 = f_inner3 + exp(- ...
				(sum(abs(eqH*(x_mod(:,l)-x_mod(:,m)) *ones(1,n_smp^4) + noise_vec).^2,1) ));
		end;
		
		f_inner = zeros(2,2);
		for num = 1:n_smp^4
			f_inner = f_inner + f_inner1(:,num)*f_inner2(num,:)/f_inner3(num);
		end;
		
		f_avg = f_avg + f_inner;
		
	end;
	MMSE = eye(2) - 1/(pi^2*2^(4*Ch.M)) * f_avg * delta^4;
	
elseif Ch.MMSE == 3
	% version 3: sum first, int second
	clear f_inner
	f_inner = zeros(2,2);
	for num = 1:n_smp^4
		for l = 1:2^(2*Ch.M)
			f_inner1 = 0;
			f_inner2 = 0;
			f_inner3 = 0;
			for k = 1:2^(2*Ch.M)
				f_inner1 = f_inner1 + x_mod(:,k)*exp(- ...
					(sum(abs(eqH*(x_mod(:,l)-x_mod(:,k)) + noise_vec(:,num)).^2,1) ));
			end;
			f_inner2 = f_inner1';
			
			for m = 1:2^(2*Ch.M)
				f_inner3 = f_inner3 + exp(- ...
					(sum(abs(eqH*(x_mod(:,l)-x_mod(:,m)) + noise_vec(:,num)).^2,1) ));
			end;
			
			f_inner = f_inner + f_inner1*f_inner2/f_inner3;
		end;
	end;
	MMSE = eye(2) - 1/(pi^2*2^(4*Ch.M)) * f_inner * delta^4;
	
elseif Ch.MMSE == 0
	MMSE = zeros(2,2);
end;