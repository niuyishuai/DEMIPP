function [diagnostics]=callipoptpoly(f,x,x0)
    % call ipopt matlab interface to solve optimization problem
    % min f(x)
    % s.t. ||x||^2<=r (r=sqrt(N))
    % with given x0
    
    if isa(f,'sdpvar') % for yalmip model
        assign(x,x0);
        N=numel(x);
        r=sqrt(N);
        output=optimize([x'*x<=r*r],f,sdpsettings('solver','ipopt','usex0',1,'verbose',0));
        diagnostics.xopt=value(x);
        diagnostics.fval=value(f);
        diagnostics.status=output.problem;
        diagnostics.iter=[];
    else
        N=numel(x0); % number of variables
        r=sqrt(N);
        grad=computegradient(f,x);
        % funcs structure
        funcs.objective = @(z)evalfcn(f,x,z);
        funcs.gradient = @(z)evalfcn(grad,x,z);
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
        diagnostics.fval=evalfcn(f,x,xopt);
        diagnostics.status=output.status;
        diagnostics.iter=output.iter;
    end
end