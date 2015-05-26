function y = logx(x)

% modified log function, which returns -inf for negtive inputs.

if length(x)>1
	error('Too many inputs. Check the logx funciton.');
end;

if x > 0
	y = log(x);
else
	y = -inf;
end;