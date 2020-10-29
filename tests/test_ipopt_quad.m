% solve
% min Jeps(x): ||x||<=r using Ipopt
%

%% build random polynomial objective function
[~,~,Q,l,s]=loadMaxCut('culberson10.txt');
N=size(Q,1); % number of variables
d=2; % degree du polynome

%% Ipopt
%close all;
% parameters
param.c=0;
param.epsilon=1e-4; % 1e-3;
param.plotflag = 1; % 1 for ploting J, only for 2D
param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function 
param.holdonflag = 0; % 1 for holding on previous plots, 0 for clean figure

%for i=1:50
y0=randonunitsphere(2*N); 
x0=y0(1:N);

[y,fval,iter,time] = proc_ipopt_quad(Q,l,s,x0,param);
round(y)
fprintf('fval: %10.5f, iter: %d, delta: %f, time: %5.2f (s).\n',fval,iter,norm(y-round(y)),time);
%end
%figure(2);view(30,30)

