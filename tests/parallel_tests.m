% Test all methods for small-scale random dataset
% min J(x): x\in {-1,1}^N

%clc;
%clear;
dirname='dataset';
filenames=dir([dirname,'/*.mat'])';
nbprobs=length(filenames);
nbmethods=6; % number of test methods
listofprobs=cell(nbprobs,1);
listofsols=cell(nbprobs,nbmethods);
listf=zeros(nbprobs,nbmethods);
listtime=listf;
listiter=listf;
listdelta=listf;
listdensity=zeros(nbprobs,1);
count=0;
coef=[]; % no use
for file=filenames
    res=regexp(file.name,'\d*\.?\d*','match');
    N=str2double(res{1});
    d=str2double(res{2});
    fprintf('Start solving problem of %d variables and of degree %d\n',N,d);
    count=count+1;
    %% load random polynomial objective function
    polytype=0; % 0: polylab (class) 1: polylab (functions), 2: yalmip, 3: syms, 4: sostools
    density=0.8; % density of polynomial
    
    % create and save large-scale data
    %[J,x,~,listofprobs{count}] = genpoly(N,d,polytype,density);
    %coef = listofprobs{count};
    %save([dirname,'/P-',num2str(N),'-',num2str(d),'.mat'],'J','x','polytype','density');
    
    % load large-scale data
    load([file.folder,'/',file.name]);
    disp(['* ',file.name]);
    J=J.simplify;
    
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
    
    %% Call houbolt scheme
    % particular parameter settings
    param.m=1;
    param.gamma=300;
    param.theta=0.8; % reduction ratio (only for variate step)
    param.tau=sqrt(2*param.m*param.epsilon);
    param.vartau = false; % true variate step, false invariate step
    param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
    
    [y,fval,iter,time] = parallel_proc(@proc_houboltscheme,J,x,N,d,coef,param);
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
    param.theta=0.8; % reduction ratio (only for variate step)
    param.tau=min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
    param.vartau = false; % true variate step, false invariate step
    param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
    
    [y,fval,iter,time] = parallel_proc(@proc_liescheme,J,x,N,d,coef,param);
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
    
    [y,fval,iter,time] = parallel_proc(@proc_rungekutta45,J,x,N,d,coef,param,tspan);
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
    [y,fval,iter,time] = parallel_proc(@proc_ipopt,J,x,N,d,coef,param);
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
    
    %% Call Quadratic reduction + Gurobi
    if 0
    [p,z]=convertbooleanpoly(J,x,0); % convert to a boolean polynomial with {0,1} variables
    z0=y0(1:N);
    [zopt,fval,iter,time] = proc_quadreduction(p,[],z0,param);
    xopt=2*zopt-1;
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
    %% Call exhaustive
    % particular parameter settings
    param.verbose=0;
    param.paralmode=1; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
        [y,fval,time]=proc_exhaustive(J,x,param);
        if param.floatingform
            fprintf('fval: %10.3f, time: %5.2f (s).\n',fval,time);
        else
            fprintf('fval: %10.3e, time: %5.2f (s).\n',fval,time);
        end
        listf(count,6)=fval;
        listiter(count,6)=iter;
        listtime(count,6)=time;
        listdelta(count,6)=0;
        listofsols{count,6}=round(y);
end

save(['parallel_results_',date,'.mat'],'listf','listiter','listtime','listdelta','listofsols','listofprobs');
