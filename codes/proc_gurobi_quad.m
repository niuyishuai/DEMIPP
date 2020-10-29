function [zopt,fval,iter,time] = proc_gurobi_quad(Q,l,s,x0);
    % gurobi for solving quadratic boolean polynomial program
    %
    
    % compute penalized objective function
    N=size(Q,1);
    
    % objective functional
    J = @(x)0.5*x'*Q*x+l'*x+s;
    
    tic;
    Q1 = 2*Q;
    id = ones(N,1);
    l1 = 2*(l - Q*id);
    s1=(0.5*id'*Q*id - l'*id + s);
    
    % solve NLP formulation (3.3) using Ipopt
    [xopt,~,iter]=callgurobimip(Q1,l1,s1,x0);
    zopt = 2*xopt-1; % convert to {-1,1}^n
    fval = J(zopt);
    time=toc;
    
    disp('Ipopt Terminated!');
end