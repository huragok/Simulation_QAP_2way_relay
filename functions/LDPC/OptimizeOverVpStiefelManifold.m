function [Vp, MI] = OptimizeOverVpStiefelManifold(Ch, lamdaX, Vp)

% Finding the right svd matrix for a fixed power allocation vector via exploiting unitary constraints.

% Reference: 
% 1. Optimization algorithms expliting unitary constraints by J. H. Manton, TSP 2002

% Author: Weiliang Zeng
% Create: 07-23-2010

if nargin == 0
	clear;
	path(path,'./func');
	
	Ch.M = 1;% 1 for BPSK; 2 for QPSK
	Ch.Ns = 2;% antenna number
	[Ch.sym_mod, Ch.sym_mod_mat] = modConfig(Ch);
	
	Ch.MI = 1;% control whether calulate MI
	Ch.MMSE = 1; % control whether and how to calulate MMSE
	
	Ch.H = sqrt(10.^(3./10))*[0.8 0.6*1i;-0.6*1i 0.8]; % channel H; for 3dB
	% Ch.G = [-1 1; 1j 1j]*[sqrt(2) 0;0 0]*[1 1j;1j 1]/2; % procoder G, the best one, maybe.
	% Ch.G = [-1 1; 1j 1j]*[sqrt(2) 0;0 0]*[1 1;1 -1]/2; % precoder G, a possible solution.
	[Ch.uh Ch.dh Ch.vh] = svd(Ch.H);
	% Ch.G = Ch.vh * eye(2) * [1 1;1 -1]/2;% inital precoder
	lamdaX = [2; 0]; % fixed value for power allocation
end;
%%
% precoding optimization
epsval = 10^-7; % stop value

maxN = 20; % the maximal iterations for steepest descent method
alpha = .3; % configuration for line search
beta = .5;

Ch.dp = diag(sqrt(lamdaX));
% Vp = [1 1j;1j 1]'/sqrt(2);% optimal right matrix. 

R = [1 0 0 0; 0 0 0 1]; % reduction matrix 
num = 1;
% fprintf('\n Optimize Right Singular Matrix:');

i = 1;
innerProdZ = 1;
%     psi = pi/6;%rand * 2 * pi;%
%     phi = pi/4;%rand * 2 * pi;%
%     Vp = [cos(psi) sin(psi)*exp(-1j*phi); -sin(psi)*exp(1j*phi) cos(psi)];
clear fVp
while innerProdZ > 10^-4
    gamma = 1;
    E = myfunMI(Ch, 0, 1, lamdaX, Vp, 1, 2);
    fVp(i) = myfunMI(Ch, 1, 0, lamdaX, Vp, 1, 2);
    MI(i) = fVp(i);
    %     fprintf('\n i = %d, fVp(i) = %f', i, fVp(i));
    %     fprintf('\n Vp=[%f+j%f, %f+j%f; %f+j%f, %f+j%f]; \n',real(Vp(1,1)),imag(Vp(1,1)),real(Vp(1,2)),imag(Vp(1,2)),...
    %         real(Vp(2,1)),imag(Vp(2,1)),real(Vp(2,2)),imag(Vp(2,2)));
    Dx = Ch.dh^2 * Ch.dp^2 * Vp * E;
    if sum(isinf(Dx))~=0
        error('Not feasible direction. Check your initianl point.')
    end;
    Gx = -Dx + Vp * Dx' * Vp;% gradient
    Z = -Gx;% descent direction
    innerProdZ = 1/2 * abs(trace(Z' * Z));
    if innerProdZ < 10^-8
        break;
    else
        while fVp(i) - myfunMI(Ch, 1, 0, lamdaX, projFunc(Vp + 2 * gamma * Z), 1, 2) >= gamma * innerProdZ
            gamma = 2 * gamma;
        end;
        while fVp(i) - myfunMI(Ch, 1, 0, lamdaX, projFunc(Vp + gamma * Z), 1, 2) < 1/2 * gamma * innerProdZ
            gamma = 1/2 * gamma;
            if gamma < 10^-6
                break
            end;
        end;
        if gamma < 10^-6
            break
        end;
        Vp = projFunc(Vp + gamma * Z);
    end;
    
    i = i + 1;
end;


%%
% plot(1:num,-MI,'o');hold on;
% xlabel('k');
% ylabel('Mutual Information');
% title('Iterative Method')


