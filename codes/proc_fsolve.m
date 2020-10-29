function [y,fval,iter,time] = proc_fsolve(J,x,y0,param)
% call fsolve for solving boolean polynomial optimization
% x\in {-1,1}^n
%
% syntax: y = proc_fsolve(J,x,y0,param)

% init parameters
N = numel(x);
epsilon = param.epsilon;
c = param.c;
options = optimoptions('fsolve','Display','none'); % using fsolve for nonlinear system

% compute gradient of the objective function J
grad = jacobian(J,x)';

tic
% Main algorithm
x0=y0(1:N);
F=@(xx)evalfcn(grad,x,xx) + c*xx + (xx.^2-1).*xx/epsilon;

[xn,~,flag,outputs]=fsolve(F,x0,options); % resoudre systeme nonlineaire pour obtenir xn+1/2
% check if the root is founded or not (the problem may have no solution)
%if flag~=1
%    disp('no solution found...');
%    return;
%end
time=toc;

fval=evalfcn(J,x,xn);
y = xn(:);
iter = outputs.iterations;
disp('Fsolve Terminated!');
end