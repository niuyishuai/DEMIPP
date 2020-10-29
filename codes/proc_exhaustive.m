function [y,fval,time] = proc_exhaustive(J,x,param)
    % Exhaustive search for solving boolean polynomial optimization
    % x\in {-1,1}^n
    %
    % syntax: y = proc_exhaustive(J,x,param)
    
    % init parameters
    N = length(x);
    verbose = param.verbose;
    paralmode = param.paralmode;
    
    tic
    % Main algorithm
    L=1:N;
    %y=ones(N,1); % case with all ones, i.e., zero -1
    %fval=evalfcn(J,x,y);
    if paralmode==0 % no parallel, sequential algorithm
        fval = inf;
        for i=0:N % case with i (i=0,...,N) -1
            fprintf('The group %d \n',i);
            allchoices=nchoosek(L,i);
            for j=1:size(allchoices,1)
                % get a vector in {-1,1}^n
                u=ones(N,1);
                u(allchoices(j,:))=-1;
                f = evalfcn(J,x,u);
                if f<fval
                    fval = f;
                    y = u;
                    if verbose
                        fprintf('%f\n',fval);
                    end
                end
            end
        end
    elseif paralmode==1 % parallel algorithm, one parfor
        fvallst=inf(N+1,1);
        ylst=zeros(N,N+1);
        %coef = param.coef;
        %d=degree(J);
        
        parfor i=1:N+1
            % recreate polynomial J for parallel process
            %x=sdpvar(N,1);
            %base=monolist(x,d);
            %J=coef'*base;
            
            fprintf('The %d-th group is started\n',i-1);
            allchoices=nchoosek(1:N,i-1);
            for j=1:size(allchoices,1)
                % get a vector in {-1,1}^n
                u=ones(N,1);
                u(allchoices(j,:))=-1;
                f = evalfcn(J,x,u);
                if f < fvallst(i)
                    fvallst(i) = f;
                    ylst(:,i) = u;
                    if verbose
                        fprintf('%f\n',fvallst(i));
                    end
                end
            end
            fprintf('The %d-th group is ended\n',i-1);
        end
        [fval,idx]=min(fvallst);
        y=ylst(:,idx);
        % elseif paralmode==2 % parallel algorithm, two parfor
        %     fvallst=inf(N,1);
        %     ylst=zeros(N,N);
        %     coef = param.coef;
        %     d=degree(J);
        %
        %     parfor i=1:N
        %         fprintf('The %d-th group is started\n',i);
        %         [fvallst(i),ylst(:,i)] = secondparfor(i,coef,N,d);
        %         fprintf('The %d-th group is finished\n',i);
        %     end
        %     [fval,idx]=min(fvallst);
        %     y=ylst(:,idx);
    end
    time=toc;
    
    disp('Exhaustive Terminated!');
end

% function [fval_ret,y_ret] = secondparfor(i,coef,N,d)
% allchoices=nchoosek(1:N,i);
% nbchoice=size(allchoices,1);
% fvallst=zeros(nbchoice,1);
% ylst=zeros(N,nbchoice);
% parfor j=1:size(allchoices,1)
%     % recreate polynomial J for parallel process
%     x=sdpvar(N,1);
%     base=monolist(x,d);
%     J=coef'*base;
%     % get a vector in {-1,1}^n
%     u=ones(N,1);
%     u(allchoices(j,:))=-1;
%     f = evalfcn(J,x,u);
%     fvallst(j)=f;
%     ylst(:,j)=u;
% end
% [fval_ret,idx]=min(fvallst);
% y_ret=ylst(:,idx);
% end
