function [xopt,fval,iter,time] = proc_ipopt(J,x,x0,param)
    % reduce a boolean polynomial function with {0,1} variables to a quadratic one,
    % then solved by GUROBI
    %
    
    % compute penalized objective function
    c=param.c;
    epsilon=param.epsilon;
    Jeps = genpenaltyfunc(J,x,c,epsilon);
    
    tic;
    % solve NLP formulation (3.3) using Ipopt
    diagnostics=callipoptpoly(Jeps,x,x0);
    xopt = diagnostics.xopt; % the computed solution (without round)
    fval = evalfcn(J,x,round(xopt)); % the value of J at round(xopt)
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
