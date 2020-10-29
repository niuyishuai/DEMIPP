% solve
% min J(x): x\in {-1,1}^N
% using Houbolt scheme
if 1
%% build random polynomial objective function
N=20; % number of variables
d=6; % degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=1; % density of polynomial
[J,x]=genpoly(N,d,polytype,density);
end
%% Main algorithm
% parameters
param.m=1; % mass
param.gamma=50; % friction factor, bigger more frictional surface
param.c=0; % for convexity

param.epsilon=1e-6; %1e-3/(4*(param.c+estimLipparam(J,x,sqrt(N)))*sqrt(N)); % 
param.explicitflag=false;
%param.tau= sqrt(2*param.m*param.epsilon);
param.tau= (3*param.gamma*param.epsilon + sqrt((3*param.gamma*param.epsilon)^2+32*param.m*param.epsilon))/4; 
param.plotflag = 1; % 1 for ploting J, only for 2D
param.surfflag = 0; % 0 no surface plot, 1 plot surface of J, 2 plot both surfaces of J and penalty function 
param.verbose = 0; % 0 scilence mode, 1 display
param.vartau = false; %true; % true variate step, false invariate step
param.theta=0.8; % reduction ratio
param.mintau = param.tau/10;%1e-2; % threshold for tau (minimal value for tau)
param.tolf = 1e-4; % tolerance for convergence of objective function
param.tolu = 1e-2; % tolerance for convergence of computed solution
param.holdonflag = 0; % 1 for holding on previous plots, 0 for clean figure

% set initial point
%y0=randonunitsphere(2*N); 
%y0=randn(2*N,1);
%y0=zeros(2*N,1);
%y0=2*rand(2*N,1)-1; 

[y,fval,iter,time] = proc_houboltscheme(J,x,y0,param);
round(y)
fprintf('fval: %10.5f, iter: %d, delta: %f, time: %5.2f (s).\n',fval,iter,norm(y-round(y)),time);    
if param.surfflag>0
    % change view angle if needed
    %figure(2);
    %view(-5,55);
end