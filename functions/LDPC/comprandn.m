function y=comprandn(M,K,deta)
% y=comprandn(M,K,deta)
% M    :  dimension
% K    :  dimension
% deta : Variance

if nargin <= 2,  
    deta = 1;
end;

if nargin <= 1,
    K = 1; 
end;

if nargin == 0,  
    M = 1; 
end;

y = (sqrt(deta/2)*randn(K,M) + 1i *sqrt(deta/2)* randn(K,M))';