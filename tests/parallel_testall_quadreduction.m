% test all methods for random problems
% min J(x): x\in {-1,1}^N

%clc;
%clear;
%setN=[10 11 12 13 14 15];
nbprocs=4;
pool=gcp('nocreate');
if isempty(pool)
    pool=parpool(nbprocs);
end

setN=[4:2:4];
setd=[2:6];
nbprobs=length(setN)*length(setd);
nbmethods=5; % number of test methods
listofprobs=cell(nbprobs,1);
listofsols=cell(nbprobs,nbmethods);
listf=zeros(nbprobs,nbmethods);
listtime=listf;
listiter=listf;
listdelta=listf;
listdensity=zeros(nbprobs,1);
if ~exist('data','dir')==0
    mkdir('data');
end
count=0;
for N=setN
    for d=setd
        fprintf('Start solving problem of %d variables and of degree %d\n',N,d);
        count=count+1;
        %% build random polynomial objective function
        polytype=0; % 0: polylab (class) 1: polylab (functions), 2: yalmip, 3: syms, 4: sostools
        density=0.8; % density of polynomial
        [J,x,~,listofprobs{count}] = genpoly(N,d,polytype,density);
        coef = listofprobs{count};
        save(['data/P-',num2str(N),'-',num2str(d),'.mat'],'J','x','polytype','density');
        
        % parameters
        param.epsilon=1e-6; % small epsilon leads to high precision
        param.theta=0.8; % reduction ratio
        param.plotflag = 0; % 1 for ploting J, only for 2D
        param.surfflag = 0; % 1 for ploting surface of J, 2 for ploting both surfaces of J and penalty function
        param.verbose = 0; % 0 scilence mode, 1 display
        param.tolf = 1e-1; % tolerance for convergence of objective function (only for Houbolt and Lie)
        
        param.c = 100;
        param.floatingform = true; % print result using floating form
        
        %% Call houbolt scheme
        % particular parameter settings
        param.m=1;
        param.gamma=100;
        param.tau=min(sqrt(2*param.m*param.epsilon),0.1);
        param.vartau = false; % true variate step, false invariate step
        param.mintau = param.tau/10;
        
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
        param.tau=min(param.epsilon/(1-param.c*param.epsilon),0.1); % tau < epsilon/(1-c*epsilon)
        param.vartau = true; % true variate step, false invariate step
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
        tspan = [0 1];
        
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
        
        %% Call Quadratic reduction + Gurobi
        [y,fval,iter,time] = parallel_proc(@proc_quadreduction,J,x,N,d,[],param);
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
        
        %% Call exhaustive
        % particular parameter settings
        param.verbose=0;
        param.paralmode=1; % used only for exhaustive method, 0: no parallel, 1: one parallel, 2: multiparallel
        %param.coef=listofprobs{count}; % used only for parallel computing (reconstruct polynomials for each worker)
        
        [y,fval,time]=proc_exhaustive(J,x,param);
        if param.floatingform
            fprintf('fval: %10.3f, time: %5.2f (s).\n',fval,time);
        else
            fprintf('fval: %10.3e, time: %5.2f (s).\n',fval,time);
        end
        listf(count,5)=fval;
        listiter(count,5)=iter;
        listtime(count,5)=time;
        listdelta(count,5)=0;
        listofsols{count,5}=round(y);
    end
end

save(['parallel_results_',date,'.mat'],'listf','listiter','listtime','listdelta','listofsols','listofprobs');

return;
delete(pool);
pool=[];

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
