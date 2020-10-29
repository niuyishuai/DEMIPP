function [y,fval,iter,time] = proc_rungekutta45_quad(Q,l,s,y0,tspan,param)
    % rungekutta45 method for solving quadratic Boolean program
    % min 0.5* x'*Q*x + l'*x + s where x in {-1,1}^n
    % syntax: y = proc_rungekutta45_quad(Q,l,s,y0,tspan,param)
    
    % init parameters
    N = size(Q,1);
    m = param.m;
    gamma = param.gamma;
    epsilon = param.epsilon;
    c = param.c;
    plotflag = param.plotflag;
    surfflag = param.surfflag;
    holdonflag = param.holdonflag;
    tolf = param.tolf;
    tolu = param.tolu;
    if plotflag
        ops=odeset('OutputFcn',@odephas2,'Events',@eventfun);
    else
        ops=odeset('Events',@eventfun);
    end
    
    % objective functional
    J = @(x)0.5*x'*Q*x+l'*x+s;
    
    % compute gradient of the objective function J
    grad = @(x)Q*x+l;
    
    % initialize plot settings
    initplot(N,plotflag,surfflag,holdonflag);
    
    % start plots
    Feps = startplot(J,N,c,epsilon,surfflag);
    
    % 时间回调函数
    function [value,isterminal,direction] = eventfun(~,y)
        U=y(1:N);
        RU=round(U); % rounding U
        Deltau=norm(U-RU);
        Deltaf=abs(J(U)-J(RU));
        if Deltaf < tolf || Deltau < tolu
            value=0; %触发时间，当其值为0的时候，时间会触发
        else
            value=1;
        end
        isterminal=1; %设为1时会，触发时间会停止求解器，设0时触发不影响工作
        direction=0; %触发方向设1时是上升触发，设-1是下降触发，设0是双向触发
    end
    
    % 使用嵌套函数参数化
    function dydt = odefcn(~,y)
        % t: time
        % y: solution
        U=y(1:N);
        P=y(N+1:end);
        dydt = [P; -(gamma*P + (U.^3-U)/epsilon + c*U + grad(U))/m];
    end
    
    tic;
    % call ode45 to solve ode
    [t,y] = ode45(@odefcn, tspan, y0, ops);
        
    % continue plots
    continueplot(J,Feps,y0,y,t,N,plotflag,surfflag);
    time=toc;
    
    y=y(end,1:N)';
    fval=evalfcn(J,[],round(y));
    iter=length(t);
    disp('Runge-Kutta45 Terminated!');
end

function initplot(N,plotflag,surfflag,holdonflag)
    % clear and initialize figures
    if plotflag
        if holdonflag == 0
            close all;
        end
        % clear surface plot
        if surfflag>0 && N==2
            figure(2);            
            hold on;
            setupfig(2);
        end
        figure(1);
        hold on
        setupfig;
    end
end

function [Feps]=startplot(J,N,c,epsilon,surfflag)
    % start plots (for first two points and the segement joining them)
    Feps = [];
    % compute penalty objective function (for 3D surface plot only)
    if surfflag==2 && N==2
        Feps = genpenaltyfunc(J,[],c,epsilon);
    end
end

function continueplot(J,Feps,y0,y,t,N,plotflag,surfflag)
    % continue ploting iteration points
    if plotflag
        % 2D figure
        plotstartpt2d(1,y0); % plot start point
        plotendpt2d(1,[y(end,1);y(end,2)]); % plot end point
    end
    
    if plotflag && surfflag>0 && N==2
        % 3D figure
        figure(2)
        % plot surface of J
        plotsurf(2,J,[],[-1.5,1.5],[-1.5,1.5]);
        % plot start point on J
        plotstartpt3d(2,[y(1,1),y(1,2),J(y(1,1:2)')]);
        if surfflag==2
            % plot surface of Fesp
            plotsurf(2,Feps,[],[-1.5,1.5],[-1.5,1.5]);
            % plot start point on Feps
            plotstartpt3d(2,[y(1,1),y(1,2),Feps(y(1,1:2)')]);
        end
        for i=1:length(t)-1
            % plot on surface of objective function J
            fval=J(y(i,1:2)');
            fval1=J(y(i+1,1:2)');
            plotsegment3d([y(i,1),y(i,2),fval],[y(i+1,1),y(i+1,2),fval1],'r-o');
            %[i,y(i,:),t(i)]
            if surfflag==2
                % plot on surface of Feps
                fvalb=Feps(y(i,1:2)');
                fvalb1=Feps(y(i+1,1:2)');
                plotsegment3d([y(i,1),y(i,2),fvalb],[y(i+1,1),y(i+1,2),fvalb1],'y-s');
            end
            drawnow limitrate
        end
        % plot end points
        plotendpt3d(2,[y(end,1),y(end,2),J(y(end,1:2)')]); % plot end point on J
        if surfflag == 2
            plotendpt3d(2,[y(end,1),y(end,2),Feps(y(end,1:2)')]); % plot end point on Feps
        end
    end
end
