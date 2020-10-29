% solve
% min J(x): x\in {-1,1}^N
% using Houbolt scheme using adapative parameters

%% build random polynomial objective function
N=10; % number of variables
d=2; % degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=1; % density of polynomial
[J,x]=genpoly(N,d,polytype,density);

%% Main algorithm
% parameters
param.c=100; % for convexity

param.epsilon=1e-6; % 1e-3;
param.theta=0.8; % reduction ratio
param.tau=sqrt(2*param.epsilon); %min(sqrt(2*param.m*param.epsilon),0.1);
param.plotflag = 1; % 1 for ploting J, only for 2D
param.surfflag = 0; % 0 no surface plot, 1 plot surface of J, 2 plot both surfaces of J and penalty function 
param.verbose = 0; % 0 scilence mode, 1 display
param.vartau = false; %true; % true variate step, false invariate step
param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
param.tolf = 1e-2; % tolerance for convergence of objective function
param.holdonflag = 1; % 1 for holding on previous plots, 0 for clean figure


%for i=1:50
%y0=randonunitsphere(2*N); 
%y0=randn(2*N,1);
%y0=zeros(2*N,1); %[0;0;0;0];
%y0=2*rand(2*N,1)-1; % [u0 p0]

[y,fval,iter,time] = proc_houboltscheme_ada(J,x,y0,param);
round(y)
fprintf('fval: %10.5f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
%disp(y);
%figure(1);view(30,30)
%end