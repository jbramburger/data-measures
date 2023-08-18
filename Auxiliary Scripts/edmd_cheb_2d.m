function K = edmd_cheb_2d(x,y,k,l)

% EDMD in a "clever" way to save memory
% The loop could be slow...

% Parameters
n = 2;
kx = nchoosek(k+n,n);
lx = nchoosek(l+n,n);

% Assemble matrices
Phi = zeros(lx,lx);
Psi = zeros(kx,lx);

% All sets of datapoints
for set = 1
    %load(sprintf('two_well_%i.mat',set),'x','y')
    nPoints = size(x,1);
    x = x(1:nPoints/2,:);
    y = y(1:nPoints/2,:);
    nPoints = nPoints / 2;
    step = 1e5;
    idx = 1:step;
    for i = 1:nPoints/step
        val_x = cheb_2d_eval(x(idx,:),l);
        val_y = cheb_2d_eval(y(idx,:),k);
        Phi = Phi + val_x.' * val_x;
        Psi = Psi + val_y.' * val_x;
        idx = idx + step;
    end
    % Now symmetrize and do it again (this assumes we know the symmetry of
    % the system, so it is cheating!)
    x(:,1) = -x(:,1);
    y(:,1) = -y(:,1);
    step = 1e5;
    idx = 1:step;
    for i = 1:nPoints/step
        val_x = cheb_2d_eval(x(idx,:),l);
        val_y = cheb_2d_eval(y(idx,:),k);
        Phi = Phi + val_x.' * val_x;
        Psi = Psi + val_y.' * val_x;
        idx = idx + step;
    end
end


% The matrix
K = Psi / Phi;
% K = sdpvar(kx,lx,'full');
% OBJ = K.*(K*Phi)-2*Psi.*K;
% options = sdpsettings('verbose', 0);
% optimize([],sum(OBJ(:)),options)
% K = value(K);