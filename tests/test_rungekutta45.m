% solve
% min J(x): x\in {-1,1}^N
% using Runge-Kutta (4,5)

%% build random polynomial objective function
N=2; % 15 number of variables
d=6; % 6 degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=0.8; % density of polynomial
[J,x]=genpoly(N,d,polytype,density);

%% Using ODE solver ode45 (Runge-Kutta 4,5)
%close all;
% parameters
param.m=1;
param.gamma=60;
param.c=100;

param.epsilon=1e-4; %1e-3/(4*(param.c+estimLipparam(J,x,sqrt(N)))*sqrt(N)); %1e-4; % 1e-3;
param.plotflag = 0; % 1 for ploting J, only for 2D
param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function 
param.holdonflag = 0; % 1 for holding on previous plots, 0 for clean figure
param.tolf = 1e-4; % tolerance for convergence of objective function
param.tolu = 1e-2; % tolerance for convergence of computed solution

tspan = [0 0.3];
%y0=randonunitsphere(2*N); 
%y0=zeros(2*N,1);
%y0=randn(2*N,1);
%y0=ones(2*N,1);
%y0=2*rand(2*N,1)-1; % [u0 p0]

[y,fval,iter,time] = proc_rungekutta45(J,x,y0,tspan,param);
round(y)
fprintf('fval: %10.5f, iter: %d, delta: %.2e, time: %5.2f (s).\n',fval,iter,norm(y-round(y)),time);    
if param.surfflag>0
    % change view angle if needed
    %figure(2);
    %view(-5,55);
end

return;
%%
% call exhaustive
param.verbose=1;
param.paralmode=0; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
%param.coef=coef; % used only for parallel computing (reconstruct polynomials for each worker)
[y1,fval1,time1]=proc_exhaustive(J,x,param);
fprintf('fval: %f, time: %5.2f (s).\n',fval1,time1);
if evalfcn(J,x,round(y)) == fval1
    fprintf('Solution quality: Global.\n');
else
    fprintf('Solution quality: Local.\n');
end    
