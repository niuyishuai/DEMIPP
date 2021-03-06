% solve
% min J(x): x\in {-1,1}^N
% using Lie's scheme
%% build random polynomial objective function
if 1
N=20; % number of variables
d=6; % 6 degree du polynome
polytype=0; % 0: polylab (class) 1: polylab (functions), 2: yalmip, 3: syms, 4: sostools
density=1; % 0.6 density of polynomial
tic
[J,x]=genpoly(N,d,polytype,density);
toc
%return
end
%% Call Lie scheme
% parameters
param.epsilon=1e-5; % 1e-3;
param.c=0; %min(10,1/param.epsilon); %100; %min(10,1/param.epsilon); % c < 1/epsilon
%r=sqrt(N);
%computec(grad,x,r);

param.explicitflag=false;
param.tau=param.epsilon/(1-param.c*param.epsilon); %param.epsilon/(1-param.c*param.epsilon); %min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
param.plotflag = 1; % 1 for ploting J, only for 2D
param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function 
param.verbose = 1; % 0 scilence mode, 1 display
param.vartau = false; %true; %false; % true variate step, false invariate step
param.theta=0.8; % reduction ratio
param.mintau = min(1e-5,param.tau/10); % threshold for tau (minimal value for tau)
param.tolf = 1e-4; % tolerance for convergence of objective function
param.tolu = 1e-2; % tolerance for convergence of computed solution
param.holdonflag = 0; % 1 for holding on previous plots, 0 for clean figure

% set initial point
%y0=randonunitsphere(2*N); 
%y0=randn(2*N,1);
%y0=zeros(2*N,1);
%y0=2*rand(2*N,1)-1;

[y,fval,iter,time] = proc_liescheme(J,x,y0,param);
%[y,fval,iter,time] = parallel_proc(@proc_liescheme,N,d,coef,param);
round(y)
fprintf('fval: %10.5f, iter: %d, delta: %f, time: %5.2f (s).\n',fval,iter,norm(y-round(y)),time);
