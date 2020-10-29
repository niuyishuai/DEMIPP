% solve
% min J(x): x\in {-1,1}^N
% using fsolve

%% build random polynomial objective function
N=2; % number of variables
d=4; % degree du polynome
x = sdpvar(N,1);
base=monolist(x,d);
nb=length(base); % length of basis
% for sparse polynomial
%density=rand;
%coef=round(20*sprand(nb,1,density)-10);
% for dense polynomial
coef=randi([-10,10],nb,1);
J=coef'*base; % polynomial objective function

%% Call fsolve
% parameters
param.epsilon=1e-4; % 1e-3;
param.c=10;

count=0;
for i=1:100
%y0=zeros(2*N,1); %[0;0;0;0];
y0=2*rand(2*N,1)-1; % [u0 p0]

[y,fval,iter,time] = proc_fsolve(J,x,y0,param);
fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
%disp(y);
if abs(norm(y)-sqrt(N))<0.1
    count=count+1;
end
end
count