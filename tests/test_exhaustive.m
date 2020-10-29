%%
% solve
% min J(x): x\in {-1,1}^N
% using Exhaustive scheme

%% build random polynomial objective function
N=10; % number of variables
d=8; % degree du polynome
polytype=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool
density=1; % density of polynomial
[J,x]=genpoly(N,d,polytype,density);

%%
param.verbose=1;
param.paralmode=1; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
param.coef=J.coef; % used only for parallel computing (reconstruct polynomials for each worker)
[y,fval,time]=proc_exhaustive(J,x,param);
fprintf('fval: %10.5f, time: %5.2f (s).\n',fval,time);
