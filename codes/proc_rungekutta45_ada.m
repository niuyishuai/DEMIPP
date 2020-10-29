function [y,fval,iter,time] = proc_rungekutta45_ada(J,x,y0,tspan,param)
    % rungekutta45_ada method for solving boolean polynomial optimization
    % x\in {-1,1}^n
    %
    % syntax: y = proc_rungekutta45_ada(J,x,y0,tspan,param)
    
    % init parameters
    N = length(x);
    m = 1;
    gamma = 3;
    iter=1;
    epsilon = param.epsilon;
    c = param.c;
    plotflag = param.plotflag;
    surfflag = param.surfflag;
    if plotflag
        ops=odeset('OutputFcn',@odephas2);
    else
        ops=odeset;
    end
    
    % compute gradient of the objective function J
    grad = computegradient(J,x);
    % compute penalized objective function if needed
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
    %Feps = J+c*x'*x/2+norm(x.^2-1)^2/(4*epsilon);
    
    % clear and initialize figures
    if plotflag
        % clear surface plot
        if surfflag>0 && N==2
            figure(2);
            if param.holdonflag == 0
                clf(2);
            end
            hold on;
        end
        figure(1);
        if param.holdonflag == 0
            clf(1);
        end
        hold on;
        box(gca,'on');
        xlabel('v_1');
        ylabel('v_2');
    end
    
    % 使用嵌套函数参数化
    function dydt = odefcn(~,y)
        % t: time
        % y: solution
        % N: number of integer variables
        % m,gamma,epsilon,c: parameters
        % grad: gradient of the objective function
        % x: decision variables of the objective function
        U=y(1:N);
        P=y(N+1:end);
        dydt = [P; -((gamma/iter)*P + (U.^3-U)/epsilon + c*U + evalfcn(grad,x,U))/m];
        iter=iter+1;
    end
    
    tic;
    % call ode45 to solve ode
    %[t,y] = ode45(@(t,y) odefcn(t,y,N,m,gamma,epsilon,c,grad,x), tspan, y0, ops);
    [t,y] = ode45(@odefcn, tspan, y0, ops);
    
    % show last 5 iterations
    %fprintf('The last 5 iterations:\n');
    %disp(y(end-5:end,:));
    %fprintf('The norm of gradient is: %e\n',norm(y(end,N+1:end)));
    
    % plot figures
    if plotflag
        plotstartpt(1,y0(1),y0(2)); % plot start point
        plotendpt(1,y(end,1),y(end,2)); % plot end point
        if surfflag>0 && N==2
            figure(2)
            plotsurf(2,J,x,[-1.5,1.5],[-1.5,1.5]);
            plotstartpt3d(2,y(1,1),y(1,2),evalfcn(J,x,y(1,1:2)')); % plot start point on J
            if surfflag==2
                plotsurf(2,Feps,x,[-1.5,1.5],[-1.5,1.5]);
                plotstartpt3d(2,y(1,1),y(1,2),evalfcn(Feps,x,y(1,1:2)')); % plot start point on Feps
            end
            for i=1:length(t)-1
                % plot on surface of objective function J
                fval=evalfcn(J,x,y(i,1:2)');
                fval1=evalfcn(J,x,y(i+1,1:2)');
                plot3([y(i,1),y(i+1,1)],[y(i,2),y(i+1,2)],[fval,fval1],'r-o','MarkerSize',2,'LineWidth',1.5);
                %[i,y(i,:),t(i)]
                
                if surfflag==2
                    % plot on surface of penalty function Feps
                    fvalb=evalfcn(Feps,x,y(i,1:2)');
                    fvalb1=evalfcn(Feps,x,y(i+1,1:2)');
                    plot3([y(i,1),y(i+1,1)],[y(i,2),y(i+1,2)],[fvalb,fvalb1],'y-s','MarkerSize',2,'LineWidth',1.5);
                end
                drawnow limitrate
            end
            plotendpt3d(2,y(end,1),y(end,2),evalfcn(J,x,y(end,1:2)')); % plot end point on J
            if surfflag == 2
                plotendpt3d(2,y(end,1),y(end,2),evalfcn(Feps,x,y(end,1:2)')); % plot end point on Feps
            end
        end
    end
    time=toc;

    y=y(end,1:N)';
    fval=evalfcn(J,x,round(y));
    iter=length(t);
    disp('Runge-Kutta45 Terminated!');
end
