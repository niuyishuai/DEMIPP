% Test MQLib MaxCut problems
% min J(x): x\in {-1,1}^N

%clc;
%clear;
%[~,~,Q,l,s]=loadMaxCut('culberson10.txt');
load('p4000_2.txt.mat');
%[J,x,W]=loadMaxCut('g000002.txt');
setN=size(Q,1);
setd=2;
nbprobs=length(setN)*length(setd);
nbmethods=6; % number of test methods
listofprobs=cell(nbprobs,1);
listofsols=cell(nbprobs,nbmethods);
listf=zeros(nbprobs,nbmethods);
listtime=listf;
listiter=listf;
listdelta=listf;
listdensity=zeros(nbprobs,1);
dirname='data_MaxCut';
if exist(dirname,'dir')==0
    mkdir(dirname);
end
count=0;
for N=setN
    for d=setd
        fprintf('Start solving MaxCut problem of %d variables and of degree %d\n',N,d);
        count=count+1;
        %% load random polynomial objective function
        % create and save small-scale data
        %[J,x,~,listofprobs{count}] = genpoly(N,d,polytype,density);
        %coef = listofprobs{count};
        %save([dirname,'/P-',num2str(N),'-',num2str(d),'.mat'],'J','x','N','d');

        % load small-scale data
        %load([dirname,'/P-',num2str(N),'-',num2str(d),'.mat']);
        %disp(['* ',dirname,'/P-',num2str(N),'-',num2str(d)]);
        
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
        %y0=randonunitsphere(2*N); % initial point
        
        %% Call houbolt scheme
        % particular parameter settings
        param.m=1;
        param.gamma=100;
        param.theta=0.8; % reduction ratio (only for variate step)
        param.tau=min(sqrt(2*param.m*param.epsilon),0.1);
        param.vartau = false; % true variate step, false invariate step
        param.mintau = param.tau/10;
        
        [y,fval,iter,time] = proc_houboltscheme_quad(Q,l,s,y0,param);
        %[y,fval,iter,time] = proc_houboltscheme(J,x,y0,param);
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
        param.epsilon=1e-8;
        param.theta=0.8; % reduction ratio (only for variate step)
        param.tau=min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
        param.vartau = false; % true variate step, false invariate step
        param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
        param.plotflag = 0;
        
        [y,fval,iter,time] = proc_liescheme_quad(Q,l,s,y0,param);
        %[y,fval,iter,time] = proc_liescheme(J,x,y0,param);
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
        tspan = [0 0.2];
        param.epsilon=1e-5;
        param.gamma=100;

        [y,fval,iter,time] = proc_rungekutta45_quad(Q,l,s,y0,tspan,param);
        %[y,fval,iter,time] = proc_rungekutta45(J,x,y0,tspan,param);
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
        param.epsilon=1e-5;
        param.gamma=100;
        
        [y,fval,iter,time] = proc_ipopt_quad(Q,l,s,x0,param);
        %[y,fval,iter,time] = proc_ipopt(J,x,x0,param);
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
        J=J.simplify;
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
        param.paralmode=0; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
        if 0
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
    end
end

%save(['Results_MaxCut_',createtimestamp,'.mat'],'listf','listiter','listtime','listdelta','listofsols');

return;