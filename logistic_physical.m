% -------------------------------------------------------------------------
% Discovering Physical Measures of the Logistic Map from Data
%
% This script is used to replicate the results of Section 5.1.1 of 
% Data-driven Discovery of Invariant Measures by Jason J. Bramburger 
% and Giovanni Fantuzzi. 
%
% The goal of the script is to use data to approximate the Lie derivative
% to then initiate a semidefinite program to identify the physical measure
% of the logistic map f(x) = 2x^2 - 1. The optimization objective is to get
% the first moment, denoted y_1, as close to the first moment,
% yOBJ, approximated from data.  
%
% Instead of the usual monomial basis, we use a Chebyshev basis to improve 
% numerical conditioning and accuracy. This requires one to use the ChebFun 
% package, which can be downloaded at: https://www.chebfun.org/download/
%
% Packages required: YALMIP, MOSEK, and ChebFun
%
% Written by J. Bramburger and G. Fantuzzi.
%
% -------------------------------------------------------------------------

% Clean workspace
clear all; 
close all; 
clc
yalmip clear
format long

% Add path to auxiliary scripts
addpath(genpath('Auxiliary Scripts')) 

%% Method Parameters
% k = maximal degree of monomials in P matrix
% ell = maximal degree of monomials in Q matrix (and number of moments)
% N = number of data points used to learn Lie derivative
% TOL = thresholding tolerance for the EDMD matrix
k = 10;
ell = 2*k;
N = 1e3;
TOL = 0;

%% Generate Synthetic Data

% Chebfun objects
T1 = chebpoly(0:k,[-1,1]);
T2 = chebpoly(0:ell,[-1,1]);

% Generate data
x = zeros(N,1);
x(1) = 0.25;
for n = 2:N
    % iterating the map 2x^2 - 1
    x(n) = 2*x(n-1)^2 - 1; 
end

% EDMD matrix
Q = T2( x(1:N-1) )'; 
P = T1( x(2:N) )'; 
K = edmd_with_thresholding(Q,P,TOL);

% Lie Derivative
L = K - eye(size(K));

%% Moment Matrices

% Moment vector
y = sdpvar(ell,1);

% Moment matrix
A = chebsdp_1d(k);
M0 = reshape(A*[1;y],[k+1,k+1]);

% Localizing moment matrix for [-1,1]
B = chebsdp_1d_locball(k);
M1 = reshape(B*[1;y],[k,k]);

%% Optimization procedure to identify physical meausres

% Approximate first moment from data
yObj = sum(x)/N;

% Solve
OBJ = ( y(1) - yObj ).^2;
sol = optimize([L*[1;y]==0, M0>=0, M1>=0], OBJ)

%% Build approximation to density with chebfun

M = A(:,1) * 2;
for j = 1:ell
    T = chebpoly(j,[-1,1]);
    M = M + A(:,j+1) .* sum(T);
end

M = full(reshape(M,[k+1,k+1]));
rho = chebfun(M\value([1;y(1:k)]),'coeffs');
xx = linspace(-1,1,10000);

% Plot results
figure(1)
plot(xx, 1./pi./sqrt(1-xx.^2),'--','Color',[1 69/255 79/255],'LineWidth',3)
hold on
plot(rho,'Color',[36/255 122/255 254/255],'LineWidth',4)
axis([-1 1 0 4])
xlabel('$x$','interpreter','latex')
legend('Exact','Discovered','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16)

%% L1 norm of cumulative distributions 
% --> used to produce data in Table 1

rhoStarX = chebfun(@(x) asin(x)/pi + 0.5,'splitting','on');
rhoX = cumsum(rho);
L1error = sum(abs(rhoStarX - rhoX))


