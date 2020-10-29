% solve
% min J(x): x\in {-1,1}^N
% using Runge-Kutta (4,5)

%% build random polynomial objective function
N=2; % 15 number of variables
d=6; % 6 degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=1; % density of polynomial
[J,x]=genpoly(N,d,polytype,density);

%% Using ODE solver ode45 (Runge-Kutta 4,5)
%close all;
% parameters
param.c=100;

param.epsilon=1e-6; % 1e-3;
param.plotflag = 1; % 1 for ploting J, only for 2D
param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function 
param.holdonflag = 1; % 1 for holding on previous plots, 0 for clean figure

%for i=1:50
tspan = [0 0.03];
y0=randonunitsphere(2*N); 
y0(N+1:end)=zeros(N,1);
%y0=[0;0;1;1];
%y0=zeros(2*N,1);
%y0=randn(2*N,1);
%y0=ones(2*N,1);
%y0=2*rand(2*N,1)-1; % [u0 p0]

[y,fval,iter,time] = proc_rungekutta45_ada(J,x,y0,tspan,param);
round(y)
fprintf('fval: %10.5f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
%figure(2);view(30,30)
%end
