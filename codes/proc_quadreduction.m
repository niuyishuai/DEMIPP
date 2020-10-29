function [y,fval,iter,time] = proc_quadreduction(J,x,x0,param)
    % reduce a boolean polynomial function with {0,1} variables to a quadratic one,
    % then solved by GUROBI
    %
    
    % init parameters
    if isa(J,'MPOLY')
        polytype=0;
        N = J.n;
    end
    
    tic;
    [~,Q,q,c]=reducemin(J,x,polytype);
    y0=zeros(length(q),1);
    y0(1:N)=x0(:);
    
    % solve Boolean quadratic problem use Gurobi
    [yopt,~,iter] = callgurobimip(Q,q,c,y0);

    % clear and initialize figures
    if param.plotflag
        % clear surface plot
        figure(1);
        if param.holdonflag == 0
            clf(1);
        end
        hold on;
        setupfig;
        
        plotsegment2d(y0,yopt,'b-o');
        plotstartpt2d(1,y0); % plot start point
        plotendpt2d(1,yopt); % plot end point
        drawnow limitrate;
    end
    y=yopt(1:N); % the computed solution (without round)
    fval=J.eval(round(y)); % the value of J at round(y)
    time=toc;
    
    disp('Quadratic reduction + Gurobi Terminated!');
end
