function [diagnostics]=callipoptpoly_quad(f,grad,x0)
    % call ipopt matlab interface to solve optimization problem
    % min f(x)
    % s.t. ||x||^2<=r (r=sqrt(N))
    % with given x0
    
    N=numel(x0); % number of variables
    r=sqrt(N);
    % funcs structure
    funcs.objective = f;
    funcs.gradient = grad;
    funcs.constraints = @(z)sum(z.^2);
    funcs.jacobian = @(z)sparse(2*z');
    funcs.jacobianstructure = @()sparse(ones(1,N));
    % opts structure
    opts.lb = -r*ones(N,1);
    opts.ub = r*ones(N,1);
    opts.cl = 0;
    opts.cu = r*r;
    opts.ipopt = ipoptset;
    opts.ipopt.hessian_approximation = 'limited-memory';
    opts.ipopt.mu_strategy           = 'adaptive';
    %opts.ipopt.tol                   = 1e-7;
    opts.ipopt.print_level           = 0;
    
    [xopt,output] = ...
        ipopt(x0,funcs,opts);
    diagnostics.xopt=xopt;
    %diagnostics.fval=f(xopt);
    diagnostics.status=output.status;
    diagnostics.iter=output.iter;
end