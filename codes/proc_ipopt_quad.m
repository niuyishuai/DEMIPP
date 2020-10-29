function [xopt,fval,iter,time] = proc_ipopt_quad(Q,l,s,x0,param)
    % opopt for solving quadratic boolean polynomial program
    %
    
    % compute penalized objective function
    c=param.c;
    epsilon=param.epsilon;
    
    % objective functional
    J = @(x)0.5*x'*Q*x+l'*x+s;
    
    % penalized function
    Jeps = @(x)J(x)+c*(x'*x)/2+(1/(4*epsilon))*(sum((x.^2-1).^2));
    
    % compute gradient of the objective function J
    grad = @(x)Q*x+l+c*x+(x.^3-x)/epsilon;
    
    tic;
    % solve NLP formulation (3.3) using Ipopt
    diagnostics=callipoptpoly_quad(Jeps,grad,x0);
    xopt = diagnostics.xopt; % the computed solution (without round)
    fval = J(round(xopt)); % the value of J at round(xopt)
    iter = diagnostics.iter;

    % clear and initialize figures
    if param.plotflag
        % clear surface plot
        figure(1);
        if param.holdonflag == 0
            clf(1);
        end
        hold on;
        setupfig;
        
        plotstartpt2d(1,x0); % plot start point
        plotendpt2d(1,xopt); % plot end point
        plotsegment2d(x0,xopt,'b-o');
        drawnow limitrate;
    end
    time=toc;
    
    disp('Ipopt Terminated!');
end
