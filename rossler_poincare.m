% -------------------------------------------------------------------------
% Discovering UPOs of the Rossler Poincare Map
%
% This script is used to replicate the results of Section 5.3.2 of 
% Data-driven Discovery of Invariant Measures by Jason J. Bramburger 
% and Giovanni Fantuzzi. 
%
% The goal of the script is to use collect data from the Poincare map x_1 =
% 0 of the Rossler system and then discover ergodic measures on that data.
% The minimizers of the ergodic measures correspond (approximately) to
% unstable periodic orbits intersected with the Poincare section which can
% then be leveraged to identify the periodic orbits embedded in the
% continuous-time attractor.
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

%% Simulate the Rossler system

% ODE45 specifications
N = 1e6;
dt = 0.005;
tspan = (0:N)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

% Simulate trajectory
x0 = [0; -15; 0]; 
rossler = @(t,x) [-x(2) - x(3); x(1) + 0.1*x(2); 0.1 + x(3)*(x(1) - 18)];
[tspan, xdat] = ode45(@(t,x) rossler(t,x),tspan,x0,options);

%% Gather Poincare section data

% Initialize data
x = [];

% Counting number of Poincare intersections
count = 1;

% Find intersections
for ind = 1:N
    if  (xdat(ind,1) < 0 && xdat(ind+1,1) >= 0)  
        x(count,:) = xdat(ind+1,2); %nth iterate
        count = count + 1;
    end
end

% Plot Poincare section data
figure(1)
plot(x(1:end),x(1:end),'Color',[0.5 0.5 0.5],'LineWidth',2)
hold on
plot(x(1:end-1),x(2:end),'k.','MarkerSize',20)
plot(-22.9179234809962,-22.9179234809962,'.','Color',[1 69/255 79/255],'MarkerSize',50) % fixed point
xlabel('$x_{2,n}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_{2,n+1}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
%title('Rossler Poincare Section Data','Interpreter','latex','FontSize',20,'FontWeight','Bold')
axis tight
set(gca,'fontsize',16)

%% Scale data into [-1,1]

% Find data extremes
xmin = min(x);
xmax = max(x);

% Scale data to [-1,1]
xs = 2*(x - xmin)/(xmax - xmin) - 1;

% Plot scaled Poincare section data
figure(2)
plot(xs(1:end-1),xs(2:end),'k.','MarkerSize',20)
xlabel('$\tilde{x}_{2,n}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$\tilde{x}_{2,n+1}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Scaled Rossler Poincare Section Data','Interpreter','latex','FontSize',20,'FontWeight','Bold')

%% Build Lie derivative

% Method Parameters
% k = maximal degree of monomials in P matrix
% ell = maximal degree of monomials in Q matrix (and number of moments)
% TOL = thresholding tolerance for the EDMD matrix
k = 20; % must be even! 
ell = 4*k;
TOL = 0;

% Chebfun objects
T1 = chebpoly(0:k,[-1,1]);
T2 = chebpoly(0:ell,[-1,1]);

% EDMD matrix
Phi = T2( xs(1:N-1) )'; 
Psi = T1( xs(2:N) )'; 
K = edmd_with_thresholding(Phi,Psi,TOL);

% Lie Derivative
L = K - eye(size(K));

%% Moment Matrices

% Moment vector
y = sdpvar(ell,1);

% Moment matrix
A = chebsdp_1d(ell/2);
M0 = reshape(A*[1;y],[ell/2+1,ell/2+1]);

% Localizing moment matrix for [-1,1]
B = chebsdp_1d_locball(ell/2);
M1 = reshape(B*[1;y],[ell/2,ell/2]);

%% Optimization procedure to identify invariant measures
% --> Loop over random objectives to identify lots of UPOs

numIters = 10;
upos = zeros(numIters,100);

for ind = 1:numIters

    OBJ = dot(10*rand(1,length(y(1:5)))-5,y(1:5)); 
    sol = optimize([L*[1;y]==0, M0>=0, M1>=0], OBJ);


    if sum(sum(isnan(value(M0)))) == 0
        mins = cheb_extractMinimizers(value(M0), 0:ell/2);
        
        % Save UPOs as | period | values of UPOs in Poincare section --
        % rest is populated by zeros
        upos(ind,1:length(mins)+1) = [length(mins), 0.5*(mins + 1)*(xmax - xmin) + xmin];
    
    end
        

end

%% Poincare Section Plots 
% --> to reproduce plots in Figure 6

% Plot Poincare section data
figure(3)
plot(x(1:end),x(1:end),'Color',[0.5 0.5 0.5],'LineWidth',2)
hold on
plot(x(1:end-2),x(3:end),'k.','MarkerSize',20)
plot(-22.9179234809962,-22.9179234809962,'.','Color',[1 69/255 79/255],'MarkerSize',50) % fixed point
plot([-18.0313858364717,-24.4815932442540],[-18.0313858364717,-24.4815932442540],'.','Color',[36/255 122/255 254/255],'MarkerSize',50) % Period 2 point
xlabel('$x_{2,n}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_{2,n+2}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
axis tight
set(gca,'fontsize',16)

figure(4)
plot(x(1:end),x(1:end),'Color',[0.5 0.5 0.5],'LineWidth',2)
hold on
plot(x(1:end-3),x(4:end),'k.','MarkerSize',20)
plot(-22.9179234809962,-22.9179234809962,'.','Color',[1 69/255 79/255],'MarkerSize',50) % fixed point
plot([-25.7885341999335,-14.0390147030195,-19.1857534188735],[-25.7885341999335,-14.0390147030195,-19.1857534188735],'.','Color',[254/255 174/255 0/255],'MarkerSize',50) % Period 3 point
plot([-16.1726475968098,-25.0761646600292,-22.0681328549753],[-16.1726475968098,-25.0761646600292,-22.0681328549753],'.','Color',[151/255 14/255 83/255],'MarkerSize',50) % Period 3 point
xlabel('$x_{2,n}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_{2,n+3}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
axis tight
set(gca,'fontsize',16)

figure(5)
plot(x(1:end),x(1:end),'Color',[0.5 0.5 0.5],'LineWidth',2)
hold on
plot(x(1:end-4),x(5:end),'k.','MarkerSize',20)
plot(-22.9179234809962,-22.9179234809962,'.','Color',[1 69/255 79/255],'MarkerSize',50) % fixed point
plot([-18.0313858364717,-24.4815932442540],[-18.0313858364717,-24.4815932442540],'.','Color',[36/255 122/255 254/255],'MarkerSize',50) % Period 2 point
plot([-16.9994950231043,-22.2068733571595,-23.1419355614561,-24.8175866399386],[-16.9994950231043,-22.2068733571595,-23.1419355614561,-24.8175866399386],'.','Color',[0 168/255 0],'MarkerSize',50) % Period 4 point
xlabel('$x_{2,n}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_{2,n+4}$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
axis tight
set(gca,'fontsize',16)

