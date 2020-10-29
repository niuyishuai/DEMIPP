% solve
% min J(x): x\in {-1,1}^N using quadratic reduction + Gurobi
%

%% build random polynomial objective function
N=5; % 15 number of variables
d=4; % 6 degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=0.8; % density of polynomial
[J,x]=genpoly(N,d,polytype,density); % J is a boolean polynomial with {-1,1} variables

%%
J=J.simplify;
[p,z]=convertbooleanpoly(J,x,0); % convert to a boolean polynomial with {0,1} variables

%% Quadratic reduction
%close all;
% parameters
param.plotflag = 0; % 1 for ploting J, only for 2D
param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function 
param.holdonflag = 1; % 1 for holding on previous plots, 0 for clean figure

for i=1:1
%z0=zeros(N,1);
%z0=2*rand(N,1)-1;
y0=randonunitsphere(2*N); 
z0=y0(1:N);

[zopt,fval,iter,time] = proc_quadreduction(p,[],z0,param);
xopt=2*zopt-1 % convert solution to {-1,1}
fprintf('fval: %10.5f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
end
%figure(2);view(30,30)

