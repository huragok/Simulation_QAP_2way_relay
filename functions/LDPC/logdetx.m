function y = logdetx(x,t)

% logdet(x) with log barrier

y = -logx(det(x)) - 1/t*( logx(2*5-trace(x)) + logx(trace(x)-2*2) );