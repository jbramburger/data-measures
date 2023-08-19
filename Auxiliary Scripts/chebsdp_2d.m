function At = chebsdp_2d(d)

% Create SDP data matrices for an SOS constraint in the chebyshev basis
% Precisely, returns coefficient matrices such that
%
% [p_00]*[p_00 p_10 ... p_0d] = A_0 p_00 + A_1 p_10 + ... + A_k p_{0,2d}
% [p_10]
% [...]
% [p_0d]
%
% The matrices A0,...A_k are stored as columns of a sparse matrix At. The
% total degree of the lower basis is d, the total degree of the upper basis
% is 2d, and the basis is products of chebyshev polynomials:
%
% p_ij(x,y) = T_i(x)*T_j(y)

% Set up 
dk = nchoosek(d+2,2);           % matrix size
A = monpowers(2,d);             % bivariate exponents up to degree d
P = monpowers(2,2*d);           % bivariate exponents up to degree 2d
row = [];
col = [];
val = [];
nB = 4*size(A,1);
hash = rand(2,1);
Phash = P*hash;
idx = (1:dk)';
for k = 1:dk
    % Exponents of product of basis vector (column) by element k of its
    % transpose: this gives the monomials of the first column of the moment
    % matrix
    B = abs([A + A(k,:); ...
             A + A(k,:).*[-1 1]; ...
             A + A(k,:).*[1 -1]; ...
             A - A(k,:)]);
    Bhash = B*hash;
    [LOA,LOCB] = ismembertol(Bhash,Phash,1e-8);
    row = [row; repmat(idx, 4, 1)];
    col = [col; LOCB];
    val = [val; 0.25*ones(nB,1)];
    idx = idx + dk;
end
At = sparse(row,col,val,dk^2,nchoosek(2*d+2,2));