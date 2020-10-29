%load ex01.mat
N=20; % number of variables
d=4; % degree du polynome

genpoly
x = sdpvar(N,1);
base=monolist(x,d);
nb=length(base); % length of basis
coef=randi([-10,10],nb,1);

J=coef'*base;
loops=100;
tic;
for i=1:loops
    x0=2*rand(N,1)-1; 
    fval=evalfcn(J,x,x0);
end
time=toc;

avgt=time/loops;
disp(avgt);
