%function [y,fval,iter,time]=test_rk45(N,d)
    % solve
    % min J(x): x\in {-1,1}^N
    % using Runge-Kutta (4,5)
    %% build random polynomial objective function
    N=6; % number of variables
    d=8; % degree du polynome
    
    tic
    polytool=0; % 0: polylab (class), 1: polymat (function), 2: yalmip, 3: syms, 4:sostool

    switch polytool
        case 0
            x = MPOLY.mpolyvars(N);
            base = MPOLY.monolist(N,d);
        case 1
            x = mpolyvars(N);
            base=mpoly_monolist(N,d);
        case 2
            x = sdpvar(N,1);
            base=monolist(x,d);
        case 3
            x = sym('x',[N,1],'real');
            base=monolist(x,d);
        case 4
            x = mpvar('x',[N,1]);
            base=monomials(x,0:d);
    end
    nb=length(base); % length of basis
    % for sparse polynomial
    %density=0.8;
    %coef=round(20*sprand(nb,1,density)-10);
    % for dense polynomial
    %coef=randi([-10,10],nb,1);
    
    switch polytool
        case 1
            J = mpoly_mtimes_leftdoublematrix(coef',base);
            %J = mpoly_simplify(J);
        otherwise
            J=coef'*base; % polynomial objective function
    end
    toc
    %return;
    %% Using ODE solver ode45 (Runge-Kutta 4,5)
    close all;
    % parameters
    param.m=1;
    param.gamma=50;
    param.c=1;
    
    param.epsilon=1e-5; % 1e-3;
    param.plotflag = 1; % 1 for ploting J, only for 2D
    param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function
    
    tspan = [0 0.5];
    %y0=[0;0;1;1];
    y0=zeros(2*N,1);
    %y0=2*rand(2*N,1)-1; % [u0 p0]
    
    [y,fval,iter,time] = proc_rungekutta45(J,x,y0,tspan,param);
    fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,int64(iter),time);
    %figure(2);view(30,30)
    
%end