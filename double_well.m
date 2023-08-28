% -------------------------------------------------------------------------
% Discovering the Invariant Measure of a Stochastic Double-Well 
%        System from Data
%
% This script is used to replicate the results of Section 6.2 of 
% Data-driven Discovery of Invariant Measures by Jason J. Bramburger 
% and Giovanni Fantuzzi. 
%
% The goal of the script is to use data from a stochastic double-well 
% system to identify the unique invariant measure associated to the system.
% The invariant measure is symmetric with two peaks centred at each well,
% but biases in the data will lead to (slightly) asymmetric measures.
%
% Computations are performed using a Chebyshev polynomial basis in 2
% dimensions.
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
k = 6; % degree of invariance constraints
ell = k + 2; % degree of EDMD dictionary (must be even)
ic = 500;  % number of initial conditions
ts = 1e4; % number of timesteps for each simulation
h = 1e-3; % timestep

%% Generate Synthetic Data
% Uses Euler-Maruyama method for simulating SDEs
% Potential is V(x,y) = (x^2 - 1) + y^2

% Noise amplitude
sigma = 0.7;

% Double-well system
f = @(x) [-2*(x(1)^2-1).*x(1); -2*x(2)];

% Initializations
m = ic*ts;
x = zeros(2,m);
y = zeros(2,m);
shift = 0;
for run = 1:ic
    x(:,shift+1) = [(-1)^mod(run,2); 0];
    y(:,shift+1) = x(:,shift+1) + h*f(x(:,shift+1)) + sqrt(h)*sigma*randn(2,1);
    for i = 1:ts-1
        x(:,shift+i+1) = y(:,shift+i);
        y(:,shift+i+1) = x(:,shift+i+1) + h*f(x(:,shift+i+1)) + sqrt(h)*sigma*randn(2,1);
    end
    shift = shift + ts;
end
x = x(:,1:m).'./2;
y = y.'./2;

% Save data for reproducability
%save('two_well_2.mat','x','y','h');

%% Plot data

% Plot a point cloud of the datapoints in x
figure(1)
plot(x(:,1),x(:,2),'k.','MarkerSize',2)
grid on

%% Approximate Lie derivative from data using EDMD

K = edmd_cheb_2d(x,y,k,ell);
Lie =  K - eye(size(K)); % Constraints are homogeneous, scale by h!
Lie = Lie(2:end,:); % remove evolution of the constant

%% Build and solve the SDP

lx = nchoosek(ell+2,2);
y = sdpvar(lx-1,1);
At = chebsdp_2d(ell/2);
ss = nchoosek(ell/2+2,2);
M = reshape(At*[1;y], ss, ss);
OBJ = -2-y(3)-y(5);
CNSTR = [M>=0, Lie*[1;y]==0];
optimize(CNSTR,OBJ)
y = value(y);

%% Extract Approximate Density
% Recall: everything is in Chebyshev basis

P1 = monpowers(2,k);
nk = nchoosek(k+2,2);
M0 = ones(nk,nk);
for j = 1:nk
    for i = 1:j
        for var = 1:2
            pow = [P1(i,var)+P1(j,var), abs(P1(i,var)-P1(j,var))];
            foo = [0, 0];
            foo(pow~=1) = ( 1 + (-1).^(pow(pow~=1)) )./(1 - (pow(pow~=1)).^2);
            M0(i,j) = M0(i,j) * (foo*[0.5; 0.5]);
        end
        M0(j,i) = M0(i,j);
    end
end
CC = zeros(k+1,k+1);
CC(sub2ind([k+1,k+1], P1(:,2)+1, P1(:,1)+1)) = M0\[1; y(1:nk-1)];
rho = chebfun2(CC,'coeffs');

%% Plot Approximate Density

h = surf(rho, 'NUMPTS', 300);

% Remove the scaling of the data
h.XData = h.XData .* 2;
h.YData = h.YData .* 2;
axis([-2 2 -2 2])
view([0 90]);
axis square
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 8;
ax.XTick = [-2 0 2];
ax.YTick = [-2 0 2];
box on
caxis([-1,3])
title(sprintf('$(k,\\ell)=(%i,%i)$',k,ell),'Interpreter','latex','FontSize',16)
xlabel('$x_1$','Interpreter','latex','FontSize',16);
ylabel('$x_2$','Interpreter','latex','FontSize',16);
set(gca,'fontsize',16)
