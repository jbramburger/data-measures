function K = edmd_with_thresholding(PsiX,PsiY,TOL)

% EDMD method with thresholding to clean small coefficients from the
% Koopman approximation

% Useful matrices
A = PsiY * PsiY.';
B = PsiX * PsiY.';
C = PsiX * PsiX.';

% Set up optimization problem
m = size(PsiX, 1);
n = size(PsiY, 1);
K = sdpvar(n, m, 'full');
options = sdpsettings('verbose', 0);
OBJ = trace(A) - 2.*trace(K*B) + trace( (K.'*K)*C );

% First solve
optimize([], OBJ, options);
REMOVE = abs( value(K) ) <= TOL;

% Iterate if necessary
CLEAN = any(REMOVE(:));
while CLEAN
    optimize(K(REMOVE)==0, OBJ, options);
    NEW_REMOVE = abs( value(K) ) <= TOL;
    CLEAN = any( NEW_REMOVE(:) - REMOVE(:) );
    REMOVE = NEW_REMOVE;
end
K = value(K);