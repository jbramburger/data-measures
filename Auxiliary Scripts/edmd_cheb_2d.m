function K = edmd_cheb_2d(x,y,k,l)

% Parameters
n = 2;
kx = nchoosek(k+n,n);
lx = nchoosek(l+n,n);

% Assemble matrices
Phi = zeros(lx,lx);
Psi = zeros(kx,lx);

% All sets of datapoints
for set = 1
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
end


% Return Koopman matrix
K = Psi / Phi;
