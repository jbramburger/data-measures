% -------------------------------------------------------------------------
% Discovering Physical Measures of the Rossler System from Data
%
% This script is used to replicate the results of Section 6.3.1 of 
% Data-driven Discovery of Invariant Measures by Jason J. Bramburger 
% and Giovanni Fantuzzi. 
%
% The goal of the script is to use data from the Rossler system to 
% approximate to the physical measure supported on the attractor.
% Optimization objectives are to fit approximations of the linear moments 
% y_100, y_010, and y_001 which are estimated from time-averages of the 
% data. 
%
% Packages required: YALMIP and MOSEK
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
% k = max degree of v_k dictionary of obserables
% ell = max degree of v_ell dictionary of obserables
k = 14;
ell = k + 1;

%% Generate Data

% Scaling factors for the variables
scale = [30; 30; 30];
shift = [0; 0; 1];

% Rossler parameters
a = 0.1;
b = 0.1;
c = 18;

% Simulate rossler (plot should be in the unit box)
dt = 0.005;
t = 0:dt:1e3;
N = length(t);
f = @(t,x)[-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];
[tt,xdat] = ode23(f,t,[0 -20 0]);
xdat = xdat./scale' - shift';

figure(1)
plot3(xdat(:,1),xdat(:,2),xdat(:,3),'k','LineWidth',1)
box on
grid on

%% Create P and Q matrices

% P matrix
pow = monpowers(3,ell);
ellx = size(pow,1); % number of nontrivial monomials in P
P = zeros(ellx,N-1);
for i = 1:ellx
   zx = xdat(1:end-1,:).^pow(i,:);
   P(i,:) = prod(zx,2);
end

% Q matrix
pow = monpowers(3,k);
kx = size(pow,1); % number of nontrivial monomials in Q
Q = zeros(kx,N-1);
for i = 1:kx
   zy = xdat(2:end,:).^pow(i,:);
   Q(i,:) = prod(zy,2)';
end

%% Create Koopman and Lie derivative matrix approximations

% Koopman approximation
K = Q*pinv(P);

% Lie approximation
L = (K - eye(size(K)))/dt;
L = L/100; % scale to improve performance 
%   -- > recall we only need elements of null space so this 
%        doesn't effect problem formulation

%% Build moment matrix

% Monomial vector, moment matrix
sdpvar x1 x2 x3
d = k;
v = monolist([x1,x2,x3],d,0);
M = v*v.';
u = monolist([x1,x2,x3],ell);
INV = L*u;

vars = getvariables([M>=0, INV==0]);
[mt,variabletype] = yalmip('monomtable');
isLinear = variabletype(vars)==0;
nonLinear = vars(~isLinear);

% Create a new sdpvar for the moments (excluding the constant) and
% substitute into the moment matrix and invariance constraint.
y = sdpvar(length(vars),1);
M = variablereplace(M,vars(isLinear),getvariables(y(isLinear)));
M = variablereplace(M,vars(~isLinear),getvariables(y(~isLinear)));
INV = variablereplace(INV,vars(isLinear),getvariables(y(isLinear)));
INV = variablereplace(INV,vars(~isLinear),getvariables(y(~isLinear)));

%% Solve optimization problem

avg1 = trapz(t,xdat(:,1))/t(end);
avg2 = trapz(t,xdat(:,2))/t(end);
avg3 = trapz(t,xdat(:,3))/t(end);
obj = ((y(1) - avg1)^2)/avg1^2 + ((y(2) - avg2)^2)/avg2^2 + ((y(3) - avg3)^2)/avg3^2;

CNSTR = [M>=0, INV==0];
opts = sdpsettings('solver','mosek');
optimize(CNSTR, obj, opts) 

%% Get density

vpows = monpowers(3,k);
x1pows = repmat(vpows(:,1),1,length(vpows(:,1))) + repmat(vpows(:,1),1,length(vpows(:,1)))';
x2pows = repmat(vpows(:,2),1,length(vpows(:,2))) + repmat(vpows(:,2),1,length(vpows(:,2)))';
x3pows = repmat(vpows(:,3),1,length(vpows(:,3))) + repmat(vpows(:,3),1,length(vpows(:,3)))';

x1int = 2*mod(x1pows+1,2)./(x1pows + 1); 
x2int = 2*mod(x2pows+1,2)./(x2pows + 1); 
x3int = 2*mod(x3pows+1,2)./(x3pows + 1);

ML = x1int.*x2int.*x3int;
rho = ML\value([1;y(1:kx-1)]);
rho = dot(rho,v);

%% Plot density

% Make density a MATLAB function
f = sdisplay(rho);
f = f{1};
f = replace(f,'*','.*');
f = replace(f,'^','.^');
f = eval(['@(x1,x2,x3)' f]);

% Evaluate density along trajectory
figure(2)
ftraj = f(xdat(:,1),xdat(:,2),xdat(:,3));
patch((xdat(1:199468,1)+shift(1))*scale(1),(xdat(1:199468,2)+shift(2))*scale(2),(xdat(1:199468,3)+shift(3))*scale(3),ftraj(1:199468),'FaceColor','none','EdgeColor','interp','LineWidth',1) % rainbow coloured trajectory!
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
zlabel('$x_3$','Interpreter','Latex')
colormap('jet')
set(gca,'Fontsize',16)
%axis tight
colorbar
grid on

%%  Density along x3 component to see rare events in the peaks (blue colour)
figure(3)
clf
patch([-10, t(1:20001) t(20001)+10],[-100; (xdat(1:20001,3)+shift(3))*scale(3); -100],[0; 0*xdat(1:20001,3); 0],[0; ftraj(1:20001); 0],'FaceColor','none','EdgeColor','interp','LineWidth',3)
xlabel('$t$','Interpreter','Latex')
ylabel('$x_3$','Interpreter','Latex')
axis([t(2) t(20001) 0 60])
set(gca,'Fontsize',16)
colormap('jet')
box on

%% Check physical measure values
% --> We compare integrating the density, the moments in the y vector from
%     the SDP, and the long-time averages from data

% clean up command window
clc

% Load in data from numerical integration of Rossler system up to t = 10^6
% to get well-approximated moments from time-averages
load Rossler_LongQuadAvgs.mat

% Linear moment predictions
fprintf(' x1 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x1*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(1)),avg1)
fprintf(' x2 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x2*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(2)),avg2)
fprintf(' x3 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x3*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(3)),avg3)

% Quadratic moment predictions
fprintf(' x1^2 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x1*x1*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(4)), mom_x1x1)
fprintf(' x1*x2 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x1*x2*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(5)), mom_x1x2)
fprintf(' x1*x3 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x1*x3*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(7)), mom_x1x3)
fprintf(' x2^2 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x2*x2*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(6)), mom_x2x2)
fprintf(' x2*x3 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x2*x3*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(8)), mom_x2x3)
fprintf(' x3^2 moment \n -------------- \n Density: %f \n SDP: %f \n Numerical: %f \n\n', value(int(int(int(x3*x3*rho,x1,-1,1),x2,-1,1),x3,-1,1)),value(y(9)), mom_x3x3)

