% Test all methods for small-scale random dataset
% min J(x): x\in {-1,1}^N

%clc;
%clear;
dirname='Beasley';
filenames=dir([dirname,'/*.mat'])';
nbprobs=length(filenames);
nbmethods=6; % number of test methods
listofprobs=cell(nbprobs,1);
listofsols=cell(nbprobs,nbmethods);
listf=zeros(nbprobs,nbmethods);
listtime=listf;
listiter=listf;
listdelta=listf;
count=0;
minN=inf;
maxN=-inf;
for file=filenames
    count=count+1;
    %% load random polynomial objective function
    
    % load data
    load([file.folder,'/',file.name]);
    disp(['* ',file.name]);
    
    if N<minN
        minN=N;
    end
    if N>maxN
        maxN=N;
    end
    
    % parameters
    param.explicitflag = true;
    param.plotflag = 0; % 1 for ploting J, only for 2D
    param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function
    param.floatingform = true; % print result using floating form
    param.verbose = 0; % 0 scilence mode, 1 display
    param.holdonflag = 0; % no holdonplot
    param.tolf = 1e-4; % tolerance for convergence of objective function
    param.tolu = 1e-2; % tolerance for convergence of computed solution
    param.c = 0;
    param.epsilon=1e-5; %estimeps(1e-4,param.c,J,x); %1e-7; % small epsilon leads to high precision
    
    % initial point
    y0=randonunitsphere(2*N); % initial point
    
    %% Call houbolt scheme
    % particular parameter settings
    param.m=1;
    param.gamma=100;
    param.tau=sqrt(2*param.m*param.epsilon);
    
    % not used params
    param.theta=0.8; % reduction ratio (only for variate step)
    param.vartau = false; % true variate step, false invariate step
    param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
    
    [y,fval,iter,time] = proc_houboltscheme_quad(Q,l,s,y0,param);
    if param.floatingform
        fprintf('fval: %10.3f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    else
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    end
    listf(count,1)=fval;
    listiter(count,1)=iter;
    listtime(count,1)=time;
    listdelta(count,1)=norm(round(y)-y);
    listofsols{count,1}=round(y);
    
    %% Call Lie scheme
    % particular parameter settings
    param.tau=min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
    
    [y,fval,iter,time] = proc_liescheme_quad(Q,l,s,y0,param);
    if param.floatingform
        fprintf('fval: %10.3f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    else
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    end
    listf(count,2)=fval;
    listiter(count,2)=iter;
    listtime(count,2)=time;
    listdelta(count,2)=norm(round(y)-y);
    listofsols{count,2}=round(y);
    
    
    %% Call HK45
    % particular parameter settings
    tspan = [0 0.3];
    
    [y,fval,iter,time] = proc_rungekutta45_quad(Q,l,s,y0,tspan,param);
    if param.floatingform
        fprintf('fval: %10.3f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    else
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    end
    listf(count,3)=fval;
    listiter(count,3)=iter;
    listtime(count,3)=time;
    listdelta(count,3)=norm(round(y)-y);
    listofsols{count,3}=round(y);
    
    %% Call Ipopt
    x0=y0(1:N);
    [y,fval,iter,time] = proc_ipopt_quad(Q,l,s,x0,param);
    if param.floatingform
        fprintf('fval: %10.3f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    else
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    end
    listf(count,4)=fval;
    listiter(count,4)=iter;
    listtime(count,4)=time;
    listdelta(count,4)=norm(round(y)-y);
    listofsols{count,4}=round(y);
    
    %% Call Gurobi
    if 0
    x0=y0(1:N);
    [xopt,fval,iter,time] = proc_gurobi_quad(Q,l,s,x0);
    if param.floatingform
        fprintf('fval: %10.3f, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    else
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
    end
    listf(count,5)=fval;
    listiter(count,5)=iter;
    listtime(count,5)=time;
    listdelta(count,5)=0;
    listofsols{count,5}=xopt;
    end
end

save(['results_MaxCut_Gyyy_',createtimestamp,'.mat'],'listf','listiter','listtime','listdelta','listofsols','minN','maxN');
