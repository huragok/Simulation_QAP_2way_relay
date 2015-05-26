function myfunOut = myfunforP(Ch, flagMI, flagMMSE, G, t, flagOpt)

% mutual information and MMSE calculation with log barrier and
% quadratic penalty.

eqH = Ch.H * G;
Ch.MMSE = flagMMSE; % control whether and how to calulate MMSE
Ch.MI = flagMI;
[I_finite, MMSE] = MI_MMSE(Ch, eqH);

if flagOpt == 1 % power allocation vector with eq. constraint.
    if flagMI == 1
        myfunOut = -I_finite - 1/t*( logx(x(1)) + logx(x(2)) );
    elseif flagMMSE == 1
        myfunOut = MMSE;
    end;

elseif flagOpt == 3 % power allocation vector with ineq. constraint.
    if flagMI == 1
        myfunOut = -I_finite - 1/t*( logx(2-trace(G*G')));
    elseif flagMMSE == 1
        myfunOut = MMSE;
    end;

elseif flagOpt == 2
    if flagMI == 1
        myfunOut = -I_finite;
    elseif flagMMSE == 1
        myfunOut = MMSE;
    end;
end;