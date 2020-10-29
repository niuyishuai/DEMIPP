function [y,fval,iter,time] = proc_houboltscheme_ada(J,x,y0,param)
    % Houbolt scheme for solving boolean polynomial optimization
    % x\in {-1,1}^n
    %
    % syntax: y = proc_houboltscheme(J,x,y0,param)
    
    % init parameters
    N = length(x);
    m = 1;
    gamma = 3;
    %gamma = param.gamma;
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
    options = optimoptions('fsolve','Display','none'); % set options of fsolve for nonlinear system
    
    % compute gradient of the objective function J
    grad = computegradient(J,x);
    % compute penalized objective function
    if surfflag==2 && N==2
        if isa(J,'struct')
            Feps = J;
            Feps = mpoly_plus(Feps, mpoly_mtimes_doublescalar(c/2, mpoly_mtimes_polymatrix(x',x)));
            w = mpoly_minus(mpoly_times(x,x),1);
            Feps = mpoly_plus(Feps, mpoly_mtimes_leftdoublematrix(1/(4*epsilon), mpoly_mtimes_polymatrix(w',w)));
            Feps=mpoly_simplify(Feps);
        else
            v = (x.^2-1);
            Feps = J+c*x'*x/2+(1/(4*epsilon))*v'*v;
        end
    end
    
    % time count for main algorithm
    tic;
    % Main algorithm
    u0=y0(1:N); % initial point for u0
    p0=y0(N+1:end); % initial point for p0
    
    % compute u1 and um1
    s=tau^2/(2*m);
    u1=(tau-s*gamma)*p0 + (1+s*(-c+(1/epsilon)))*u0 - (s/epsilon)*u0.^3 - s*evalfcn(grad,x,u0);
    um1 = u1-2*tau*p0;
    
    % start loop
    Un=u1;
    Unm1=u0;
    Unm2=um1;
    
    % clear and initialize figures
    if plotflag
        figure(2);
        if param.holdonflag == 0
            clf(2);
        end
        hold on
        setupfig;
        % plot u0
        plotstartpt(2,u0(1),u0(2));
        % joint u0 and u1
        plot([u0(1),u1(1)],[u0(2),u1(2)],'b-o','Color',[0 0.447 0.741]);
        % clear surface plot
        if surfflag>0 && N==2
            figure(1);
            if param.holdonflag == 0
                clf(1);
            end
            setupfig;
            hold on;
        end
    end
    
    % surface plots only for N=2
    if plotflag && surfflag>0 && N==2
        plotsurf(1,J,x,[-1.5,1.5],[-1.5,1.5]);
        plotstartpt3d(1,u0(1),u0(2),evalfcn(J,x,u0));
        if surfflag==2
            plotsurf(1,Feps,x,[-1.5,1.5],[-1.5,1.5]);
            plotstartpt3d(1,u0(1),u0(2),evalfcn(Feps,x,Un));
        end
    end
    
    Delta=inf;
    iter = 3;
    while(Delta>tolf)
        iter=iter+1;
        if vartau
            if tau > mintau
                tau = tau*theta; % update tau
            end
        end
        gamma=3/(iter+1); % using a(t)=3/t, m=1
        F=@(xx)(m/tau^2)*(2*xx - 5*Un + 4*Unm1 - Unm2) + (gamma/(2*tau))*(3*xx-4*Un + Unm1) + ...
            (1/epsilon)*(xx.^3 - xx) + c*(2*Un - Unm1) + evalfcn(grad,x,2*Un-Unm1);
        [Unp1,fval,flag]=fsolve(F,Un,options); % resoudre systeme nonlineaire pour obtenir Un+1
        % check if the root is correctly founded or not
        %if flag~=1
        %    disp('no solution found...');
        %    return;
        %end
        % compute delta
        Delta = norm(Unp1-Un,inf); %abs(fval1-fval); %norm(Unp1-Un);
        if verbose
            fprintf('%5d \t %.5e \t %.5e\n',iter,fval1, Delta);
        end
        if plotflag
            % plot phase2 figure of the first two components
            figure(2)
            plot([Un(1),Unp1(1)],[Un(2),Unp1(2)],'b-o','Color',[0 0.447 0.741]);
            drawnow limitrate;
        end
        if surfflag>0 && N==2
            figure(1)
            % plot on surface of objective function J
            fval=evalfcn(J,x,Un);
            fval1=evalfcn(J,x,Unp1);
            plot3([Un(1),Unp1(1)],[Un(2),Unp1(2)],[fval,fval1],'r-o','MarkerSize',2,'LineWidth',1.5);
            if surfflag==2
                % plot on surface of penalty function Feps
                fvalb=evalfcn(Feps,x,Un);
                fvalb1=evalfcn(Feps,x,Unp1);
                plot3([Un(1),Unp1(1)],[Un(2),Unp1(2)],[fvalb,fvalb1],'y-s','MarkerSize',2,'LineWidth',1.5);
            end
            drawnow limitrate;
        end
        Unm2=Unm1;
        Unm1=Un;
        Un=Unp1;
    end
    % plot end point
    if plotflag
        plotendpt(2,Un(1),Un(2));
        if surfflag >0 && N==2
            plotendpt3d(1,Un(1),Un(2),evalfcn(J,x,Un));
            if surfflag == 2
                plotendpt3d(1,Un(1),Un(2),evalfcn(Feps,x,Un));
            end
        end
    end
    time=toc;
    
    %fval=evalfcn(J,x,Un);
    %y = Un(:);
    y = Un(:);
    fval=evalfcn(J,x,round(y));
    disp('Adaptive Houbolt Scheme Terminated!');
end