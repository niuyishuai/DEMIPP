% test all methods for random problems
% min J(x): x\in {-1,1}^N

%clc;
%clear;
%setN=[10 11 12 13 14 15];
pool=gcp('nocreate');
if isempty(pool)
    pool=parpool(4);
end
setN=2:2:16;
setd=[4 6];
nbprobs=length(setN)*length(setd);
nbmethods=4; % number of test methods
listofprobs=cell(nbprobs,1);
listofsols=cell(nbprobs,nbmethods);
listf=zeros(nbprobs,nbmethods);
listtime=listf;
listiter=listf;
listdelta=listf;
listdensity=zeros(nbprobs,1);
count=0;
for N=setN
    for d=setd
        fprintf('Start solving problem of %d variables and of degree %d\n',N,d);
        count=count+1;
        %% build random polynomial objective function
        polytype=4; % 0: polylab (class) 1: polylab (functions), 2: yalmip, 3: syms, 4: sostools
        density=1; % density of polynomial
        [J,x,~,listofprobs{count}] = genpoly(N,d,polytype,density);
        coef = listofprobs{count};
        
        % parameters
        param.epsilon=1e-5; % 1e-3;
        param.theta=0.8; % reduction ratio
        param.plotflag = 0; % 1 for ploting J, only for 2D
        param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function
        param.verbose = 0; % 0 scilence mode, 1 display
        param.tolf = 1e-4; % tolerance for convergence of objective function
        
        param.c=100;
        y0=zeros(2*N,1);%[0;0;0;0];
        %y0=2*rand(2*N,1)-1; % [u0 p0]
        
        %% Call Lie scheme
        % particular parameter settings
        param.tau=min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
        param.vartau = false; % true variate step, false invariate step
        param.mintau = param.tau/10; % threshold for tau (minimal value for tau)
        
        [y,fval,iter,time] = parallel_proc(@proc_liescheme,J,x,N,d,coef,param);
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
        listf(count,1)=fval;
        listiter(count,1)=iter;
        listtime(count,1)=time;
        listdelta(count,1)=norm(round(y)-y);
        listofsols{count,1}=y;
        
        %% Call houbolt scheme
        % particular parameter settings
        param.m=1;
        param.gamma=100;
        param.tau=min(sqrt(2*param.m*param.epsilon),0.1);
        
        [y,fval,iter,time] = parallel_proc(@proc_houboltscheme,J,x,N,d,coef,param);
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
        listf(count,2)=fval;
        listiter(count,2)=iter;
        listtime(count,2)=time;
        listdelta(count,2)=norm(round(y)-y);
        listofsols{count,2}=y;
        
        %% Call HK45
        % particular parameter settings
        tspan = [0 1];
        
        [y,fval,iter,time] = parallel_proc(@proc_rungekutta45,J,x,N,d,coef,param,tspan);
        fprintf('fval: %10.3e, iter: %d, time: %5.2f (s).\n',fval,iter,time);
        listf(count,3)=fval;
        listiter(count,3)=iter;
        listtime(count,3)=time;
        listdelta(count,3)=norm(round(y(1:N))-y(1:N));
        listofsols{count,3}=y;
        
        %% Call exhaustive
        % particular parameter settings
        param.verbose=0;
        param.paralmode=1; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
        param.coef=listofprobs{count}; % used only for parallel computing (reconstruct polynomials for each worker)
        
        [y,fval,time]=proc_exhaustive(J,x,param);
        fprintf('fval: %f, time: %5.2f (s).\n',fval,time);
        listf(count,4)=fval;
        listiter(count,4)=iter;
        listtime(count,4)=time;
        listdelta(count,4)=0;
        listofsols{count,4}=y;
    end
end

save('parallel_results_01.mat','listf','listiter','listtime','listdelta','listofsols','listofprobs');

delete(pool);

return
%%
count=0;
figure;
hold on;
for N=setN
    for d=setd
        count=count+1;
        %phi(count)=param.epsilon*N*d;
        phi = param.epsilon*sqrt(N);
        scatter(listdelta(count,1),phi,'ro');
        scatter(listdelta(count,2),phi,'bs');
        scatter(listdelta(count,3),phi,'mo');
        xlabel('\delta');
        ylabel('Nd\epsilon');
    end
end
