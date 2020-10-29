function [y,fval,iter,time] = proc_liescheme(J,x,y0,param)
    % Lie scheme for solving boolean polynomial optimization
    % x\in {-1,1}^n
    %
    % syntax: y = proc_liescheme(J,x,y0,param)
    
    % init parameters
    N = length(x);
    epsilon = param.epsilon;
    c = param.c;
    tau = param.tau;
    theta = param.theta; % reducation factor
    plotflag = param.plotflag;
    surfflag = param.surfflag;
    verbose = param.verbose;
    vartau = param.vartau;
    mintau = param.mintau;
    tolf = param.tolf;
    tolu = param.tolu;
    holdonflag = param.holdonflag;
    explicitflag = param.explicitflag;
    options = optimoptions('fsolve','Display','none'); % set options of fsolve for nonlinear system
    p=epsilon/tau + epsilon*c - 1; % for Cardano's formula
    
    % compute gradient of the objective function J
    grad = computegradient(J,x);
    
    
    tic;
    % Main algorithm
    % initialize plot settings
    initplot(y0,N,plotflag,surfflag,holdonflag);
    
    % start plots
    [Feps,fval1] = startplot(J,x,y0,N,c,epsilon,plotflag,surfflag);
    
    % start loop
    xn=y0(1:N);
    Deltaf=inf;
    Deltau=inf;
    iter = 0;
    while(Deltaf>tolf && Deltau > tolu)
        iter=iter+1;
        if vartau
            if tau > mintau
                tau = tau*theta; % update tau
            end
        end
        % step 1
        F=@(xx)xx+tau*evalfcn(grad,x,xx)-xn;
        [xn2,~,flag]=fsolve(F,xn,options); % resoudre systeme nonlineaire pour obtenir xn+1/2
        % check if the root is founded or not (the problem may have no solution)
        %if flag~=1
        %    disp('no solution found...');
        %    return;
        %end
        % step 2
        if explicitflag == false
            % using fsolve
            phi=@(xx)(1-tau/epsilon + tau*c)*xx + tau/epsilon*xx.^3 - xn2;
            xn3=fsolve(phi,xn2,options);
        else
            % using Cardano's formula
            q=-xn2*epsilon/tau;
            qo2 = -q/2;
            D = sqrt(qo2.^2 + (p/3)^3);
            xn3 = nthroot(qo2 + D,3) + nthroot(qo2 - D,3);
        end
        % compute delta
        Deltau = norm(xn3-xn);
        fval=fval1;
        fval1=evalfcn(J,x,xn3);
        Deltaf = abs(fval1-fval); % absolute error
        %Deltaf = abs(fval1-fval)/(1+abs(fval)); % relative error
        if verbose
            fprintf('%5d \t %.5e \t %.5e \t %.5e\n',iter,fval1,Deltaf,Deltau);
        end
        
        % continue plots
        [fvalb1]=continueplot(Feps,x,xn,xn3,fval,fval1,N,plotflag,surfflag);
        
        % next iteration
        xn=xn3;
    end
    % end plot
    endplot(xn3,fval1,fvalb1,N,plotflag,surfflag);
    
    time=toc;
    
    y = xn(:);
    fval=evalfcn(J,x,round(y));
    disp('Lie Scheme Terminated!');
end

function initplot(y0,N,plotflag,surfflag,holdonflag)
    % clear and initialize figures
    if plotflag
        figure(1);
        if holdonflag == 0
            clf(1);
        end
        hold on
        setupfig;
        % plot y0
        plotstartpt2d(1,y0);
        % clear surface plot
        if surfflag>0 && N==2
            figure(2);
            if holdonflag == 0
                clf(2);
            end
            hold on;
            setupfig(2);
        end
    end
end

function [Feps,fval1]=startplot(J,x,y0,N,c,epsilon,plotflag,surfflag)
    % start plots (for first two points and the segement joining them)
    Feps = [];
    fval1=evalfcn(J,x,y0(1:N));
    % compute penalty objective function (for 3D surface plot only)
    if surfflag==2 && N==2
        Feps = genpenaltyfunc(J,x,c,epsilon);
    end
    % plot surfaces
    if plotflag && surfflag>0 && N==2
        plotsurf(2,J,x,[-1.5,1.5],[-1.5,1.5]);
        plotstartpt3d(2,[y0(1),y0(2),fval1]);
        if surfflag==2
            plotsurf(2,Feps,x,[-1.5,1.5],[-1.5,1.5]);
            plotstartpt3d(2,[y0(1),y0(2),evalfcn(Feps,x,y0(1:N))]);
        end
        drawnow limitrate;
    end
end

function [fvalb1]=continueplot(Feps,x,xn,xn3,fval,fval1,N,plotflag,surfflag)
    % continue ploting iteration points
    fvalb1=inf;
    if plotflag
        % 2D figure
        % plot phase2 figure of the first two components
        figure(1)
        plotsegment2d(xn,xn3,'b-o');
        drawnow limitrate;
    end
    % surface plots only for N=2
    if surfflag>0 && N==2
        figure(2);
        %fval=evalfcn(J,x,y0(1:2));
        plotsegment3d([xn;fval],[xn3;fval1],'r-o');
        if surfflag==2
            fvalb=evalfcn(Feps,x,xn);
            fvalb1=evalfcn(Feps,x,xn3);
            plotsegment3d([xn;fvalb],[xn3;fvalb1],'y-s');
        end
        drawnow limitrate
    end
end

function endplot(xn3,fval1,fvalb1,N,plotflag,surfflag)
    % plot end point
    if plotflag
        plotendpt2d(1,xn3);
        if surfflag >0 && N==2
            plotendpt3d(2,[xn3;fval1]); % plot end point on J
            if surfflag == 2
                plotendpt3d(2,[xn3;fvalb1]); % plot end point on Feps
            end
        end
        drawnow limitrate;
    end
end
