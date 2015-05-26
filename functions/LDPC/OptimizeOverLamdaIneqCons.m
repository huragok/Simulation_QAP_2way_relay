function [x, funMI] = OptimizeOverLamdaIneqCons(Ch, x, Vp)

% Finding the power allocation via optimizing
% Barrier method for inequality constraint 

% Reference: 
% 1. Convex Optimization by S. Boyd
% 2. Numerical optimization by J Nocedal, SJ Wright.

% Author: Weiliang Zeng
% Create: 07-22-2010

if (0)
    clc;
    clear;
    path(path,'./func');

    Ch.M = 2;% 1 for BPSK; 2 for QPSK
    Ch.Ns = 2;% antenna number
    [Ch.sym_mod, Ch.sym_mod_mat] = modConfig(Ch);

    Ch.MI = 1;% control whether calulate MI
    Ch.MMSE = 1; % control whether and how to calulate MMSE
    
    Ch.H = sqrt(10.^(3./10))*[0.8 0.6*1i;-0.6*1i 0.8]; % channel H; for 3dB
    [Ch.uh Ch.dh Ch.vh] = svd(Ch.H);
	x = [.5; .5]; % inital value for power allocation.
	psi = pi/6;%
	phi = pi/4;%
	Vp = [cos(psi) sin(psi)*exp(-1j*phi); -sin(psi)*exp(1j*phi) cos(psi)];
    Vp = [1 1j;1j 1]/sqrt(2);% maybe the best one.
end;

% precoding optimization

maxN = 12; % the maximal iterations for steepest descent method
alpha = .3; % configuration for line search
beta = .5;

epsT = 10^-6; % stop criterion for barrier method
mu = 20; % for central barrier method
t = 20; % for barrier method

R = [1 0 0 0; 0 0 0 1]; % reduction matrix 
i = 1;
fprintf('\n Optimize Power Allocation:');

while 2/t > epsT % stop criterion for central path
    n = 1;
    while n < maxN % solve for a special t
        % cal G based on a given x (W)
		
        E = myfunMI(Ch, 0, 1, x, Vp, t, 3);
        Dx = -R*reshape(Ch.dh^2 * Vp * E * Vp', 4, 1) - 1/t * (1./x - 1/ (2-sum(x)));  % descent direction;
        if sum(isinf(Dx))~=0
            error('Not feasible direction. Check your initianl point.')
        end;
        deta_x = -real(Dx);
        % backtracking line search pp. 464
        ls_t = .05;
        funMI(i) = myfunMI(Ch, 1, 0, x, Vp, t, 3);
        if isinf(funMI(i))
            error('Not feasible maybe. Check your initinal point.');
        end;
        % 		fprintf('\n\n lamdaX = [%f; %f]; ',x(1), x(2));
        % 		fprintf('\n n = %d, i = %d, funMI(i) = %f', n, i, funMI(i));
  
        %         funMI_nxt = myfunMI(Ch, 1, 0, x+ls_t*deta_x, Vp, t, 3);
        %         while myfunMI(Ch, 1, 0, x+ls_t*deta_x, Vp, t, 3) > funMI(i)...
        %                 + alpha*ls_t*Dx'*deta_x
        %             if ls_t < 10^-4
        %                 break;
        %             else
        %                 ls_t = beta * ls_t;
        %             end;
        %         end;
        
        while funMI(i) - myfunMI(Ch, 1, 0, x+ 2 * ls_t*deta_x, Vp, t, 3) >= 0 %ls_t * (deta_x') * deta_x
            ls_t = 2 * ls_t;
        end;
        while funMI(i) - myfunMI(Ch, 1, 0, x + ls_t * deta_x, Vp, t, 3) < 0 %1/2 * ls_t * (deta_x') * deta_x
            ls_t = 1/2 * ls_t;
            if ls_t < 10^-6
                break
            end;
        end;
        
        if ls_t < 10^-6
            break;
        end;
        x = x + ls_t*deta_x;
        
        n = n + 1;
        i = i + 1;
        if min(x) < .08
            x = [1.99;0.001];
            funMI(i) = myfunMI(Ch, 1, 0, x, Vp, 10e6, 3);
            break;
        end;
    end;
    if min(x) < .08
        break;
    end;
    t = mu * t;

end;
fprintf('\n\n lamdaX = [%f; %f]; ',x(1), x(2));
fprintf('\n n = %d, i = %d, funMI(end) = %f', n, i, funMI(end));
myfunMI(Ch, 1, 0, [1.99;0.001], Vp, 10e6, 3),

% 
% figure(1);
% plot(1:length(funMI),-real(funMI));
% xlabel('k');
% ylabel('Mutual Information');
% title('Iterative Method')


