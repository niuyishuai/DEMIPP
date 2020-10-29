function [y,fval,iter,time] = proc_houboltscheme_quad(Q,l,s,y0,param)
    % Houbolt scheme for solving boolean quadratic optimization
    % x\in {-1,1}^n
    %
    % syntax: y = proc_houboltscheme_quad(Q,l,s,y0,param)
    
    % init parameters
    N = size(Q,1);
    m = param.m;
    gamma = param.gamma;
    epsilon = param.epsilon;
    c = param.c;
    plotflag = param.plotflag;
    surfflag = param.surfflag;
    tau = param.tau;
    vartau = param.vartau;
    mintau = param.mintau;
    theta = param.theta;
    verbose = param.verbose;
    tolf = param.tolf; % tolerence for objective value
    tolu = param.tolu; % tolerence for optimal solution
    holdonflag = param.holdonflag;
    p=(2*m/tau + 3*gamma/2)*epsilon/tau - 1; % for Cardano's formula
    
    % objective functional
    J = @(x)0.5*x'*Q*x+l'*x+s;
    
    % compute gradient of the objective function J
    grad = @(x)Q*x+l;
    
    % time count for main algorithm
    tic;
    % Main algorithm
    u0=y0(1:N); % initial point for u0
    p0=y0(N+1:end); % initial point for p0
    
    % compute u1 and um1
    s=tau^2/(2*m);
    u1=(tau-s*gamma)*p0 + (1+s*(-c+(1/epsilon)))*u0 - (s/epsilon)*u0.^3 - s*evalfcn(grad,[],u0);
    um1 = u1-2*tau*p0;
    
    % initialize plot settings
    initplot(u0,u1,N,plotflag,surfflag,holdonflag)
    
    % start plots
    [Feps,fval1] = startplot(J,[],u0,u1,N,c,epsilon,plotflag,surfflag);
    
    Deltaf=inf;
    Deltau=inf;
    iter = 3;
    % start loop
    Un=u1;
    Unm1=u0;
    Unm2=um1;
    while(Deltaf>tolf && Deltau > tolu)
        iter=iter+1;
        if vartau
            if tau > mintau
                tau = tau*theta; % update tau
            end
        end
        % using Cardano's formula
        C=(m/tau^2)*(- 5*Un + 4*Unm1 - Unm2) + (gamma/(2*tau))*(-4*Un + Unm1) + c*(2*Un - Unm1);
        P=2*Un-Unm1;
        q=epsilon*(C + grad(P));
        qo2=-q/2;
        D=sqrt(qo2.^2 + (p/3)^3);
        Unp1= nthroot(qo2 + D,3) + nthroot(qo2 - D,3);
        % compute delta
        Deltau = norm(Unp1-Un);
        fval=fval1; %evalfcn(J,x,Un);
        fval1=J(Unp1);
        Deltaf = abs(fval1-fval);
        if verbose
            fprintf('%5d \t %.5e \t %.5e \t %.5e\n',iter, fval1, Deltaf, Deltau);
        end
        % continue plots
        [fvalb1]=continueplot(Feps,[],Un,Unp1,fval,fval1,N,plotflag,surfflag);
        
        % next iteration
        Unm2=Unm1;
        Unm1=Un;
        Un=Unp1;
    end
    % end plot
    endplot(Unp1,fval1,fvalb1,N,plotflag,surfflag)
    
    time=toc;
    
    y = Unp1(:);
    fval=J(round(y)); % compute objective value at rounding point
    disp('Houbolt Scheme Terminated!');
end

function initplot(u0,u1,N,plotflag,surfflag,holdonflag)
    % clear and initialize figures
    if plotflag
        figure(1);
        if holdonflag == 0
            clf(1);
        end
        hold on
        setupfig;
        % plot u0
        plotstartpt2d(1,u0);
        % joint u0 and u1
        plotsegment2d(u0,u1,'b-o');
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

function [Feps,fval1]=startplot(J,x,u0,u1,N,c,epsilon,plotflag,surfflag)
    % start plots (for first two points and the segement joining them)
    Feps = [];
    fval1=evalfcn(J,x,u1);
    % compute penalty objective function (for 3D surface plot only)
    if surfflag==2 && N==2
        Feps = genpenaltyfunc(J,x,c,epsilon);
    end
    % surface plots only for N=2
    if plotflag && surfflag>0 && N==2
        plotsurf(2,J,x,[-1.5,1.5],[-1.5,1.5]);
        plotstartpt3d(2,[u0;evalfcn(J,x,u0)]);
        fval=evalfcn(J,x,u0);
        plotsegment3d([u0;fval],[u1;fval1],'r-o');
        if surfflag==2
            plotsurf(2,Feps,x,[-1.5,1.5],[-1.5,1.5]);
            plotstartpt3d(2,[u0;evalfcn(Feps,x,u0)]);
            fvalb=evalfcn(Feps,x,u0);
            fvalb1=evalfcn(Feps,x,u1);
            plotsegment3d([u0;fvalb],[u1;fvalb1],'y-s');
        end
    end
end

function [fvalb1]=continueplot(Feps,x,Un,Unp1,fval,fval1,N,plotflag,surfflag)
    % continue ploting iteration points
    fvalb1=inf;
    if plotflag
        % plot phase2 figure of the first two components
        figure(1)
        plotsegment2d(Un,Unp1,'b-o');
        drawnow limitrate;
    end
    if surfflag>0 && N==2
        figure(2)
        % plot on surface of objective function J
        plotsegment3d([Un;fval],[Unp1;fval1],'r-o');
        if surfflag==2
            % plot on surface of penalty function Feps
            fvalb=evalfcn(Feps,x,Un);
            fvalb1=evalfcn(Feps,x,Unp1);
            plotsegment3d([Un;fvalb],[Unp1;fvalb1],'y-s');
        end
        drawnow limitrate;
    end
end

function endplot(Unp1,fval1,fvalb1,N,plotflag,surfflag)
    % plot end point
    if plotflag
        plotendpt2d(1,Unp1);
        if surfflag >0 && N==2
            plotendpt3d(2,[Unp1;fval1]);
            if surfflag == 2
                plotendpt3d(2,[Unp1;fvalb1]);
            end
        end
    end
end