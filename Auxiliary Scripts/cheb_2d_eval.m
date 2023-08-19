function [VAL,DEG] = cheb_2d_eval(X,d)

n = 2;
DEG = monpowers(n,d);
nDEG = size(DEG,1);
nX = size(X,1);
VAL = ones(nX,nDEG);
for i = 1:nDEG
    for j = 1:n
        C = zeros(DEG(i,j)+1, 1); C(1) = 1;
        POLYVAL = chebpolyval(C, X(:,j));
        VAL(:,i) = VAL(:,i) .* POLYVAL;
    end
end