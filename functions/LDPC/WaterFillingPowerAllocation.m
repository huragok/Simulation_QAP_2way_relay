function p = WaterFillingPowerAllocation(Ch)

% the specific alg. for water filling
% Ref.: Convex Opt. by S. Boyd p245

lambdaH = diag(Ch.dh).^2;

p = zeros(length(lambdaH),1);
mu_l = 0*ones(2,1);
mu_h = 2000*ones(2,1);
mu = (mu_l + mu_h)/2;
i = 1;
while abs(sum(p)-2)>10e-6
    
    p = (mu - 1./lambdaH > 0).*(mu - 1./lambdaH);
    
    if sum(p)<2
        mu_l = mu;
    elseif sum(p)>= 2
        mu_h = mu;
    end;
    mu = (mu_l + mu_h)/2;
    i = i + 1;
	if i > 500
		error('Maybe need to raise up the upper limit')
		break;
	end;
end;