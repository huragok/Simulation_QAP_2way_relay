function projY = projFunc(X)

% projection function

[uX dX vX] = svd(X);
projY = uX * vX';