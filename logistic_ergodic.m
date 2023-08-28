% -------------------------------------------------------------------------
% Discovering Ergodic Measures of the Logistic Map from Data
%
% This script is used to replicate the results of Section 6.1.2 of 
% Data-driven Discovery of Invariant Measures by Jason J. Bramburger 
% and Giovanni Fantuzzi. 
%
% The goal of the script is to use data to approximate the Lie derivative
% to then initiate a semidefinite program to identify ergodic measures
% of the logistic map f(x) = 2x^2 - 1 supported on the unstable fixed 
% points. The optimization objective can be any linear combination of 
% moments.  
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
k = 5;
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

%% Optimization procedure to identify ergodic meausres

OBJ = y(1); % Optimization objective - linear for ergodic measures
sol = optimize([L*[1;y]==0, M0>=0, M1>=0], OBJ );

%% Build approximation to density with chebfun and plot it

% Get density
M = A(:,1) * 2;
for j = 1:ell
    T = chebpoly(j,[-1,1]);
    M = M + A(:,j+1) .* sum(T);
end
M = full(reshape(M,[k+1,k+1]));
rho = chebfun(M\value([1;y(1:k)]),'coeffs');

% Space points
xx = linspace(-1,1,10000);

% Plot result
figure(1)
plot(xx,rho(xx),'Color',[1 69/255 79/255],'LineWidth',2);

%% Print minizers
% --> uses cheb_extractMinimizers.m which is a modified extractMinimizers 
%      for the Chebyshev basis 

xopt = cheb_extractMinimizers(value(M0), 0:ell/2)
hold on
plot(xopt,rho(xopt),'ko','MarkerSize',8,'MarkerFaceColor','k')

%% Reproduce Figure 2 in paper

% Load in ergoduc measure data with objective = y(1) (first moment)
load LogisticErgodic.mat

figure(2)
plot(xx,logErg5,'--','Color',[0 168/255 0],'LineWidth',3)
hold on 
plot(xx,logErg10,'-.','Color',[1 69/255 79/255],'LineWidth',3)
plot(xx,logErg20,'Color',[36/255 122/255 254/255],'LineWidth',3)
xlabel('$x$','interpreter','latex')
legend('$k = 5$','$k = 10$','$k = 20$','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16)
axis tight

