function [yopt,fopt,avgiter,totaltime] = parallel_proc(fhandle,J,x,N,d,coef,param,tspan)
    tic;
    spmd
        %x = sdpvar(N,1);
        %base=monolist(x,d);
        %J=coef'*base;
        %if labindex == 1
        %    x0=zeros(2*N,1);
        %else
        %    x0=2*rand(2*N,1)-1;
        %end
        % the initial point are generated randomly for each parallel run.
        x0=randonunitsphere(2*N);
        if isequal(fhandle, @proc_rungekutta45)
            [y,fval,iter] = fhandle(J,x,x0,tspan,param);
        elseif isequal(fhandle, @proc_quadreduction)
            [p,z]=convertbooleanpoly(J,x,0); % convert to a boolean polynomial with {0,1} variables
            z0=x0(1:N);
            [zopt,fval,iter] = fhandle(p,z,z0,param);
            y=2*zopt-1; % convert back to {-1,1}
        elseif isequal(fhandle, @proc_ipopt)
            y0=x0(1:N);
            [y,fval,iter] = fhandle(J,x,y0,param);
        else
            [y,fval,iter] = fhandle(J,x,x0,param);
        end
    end
    totaltime=toc;
    
    nbres=length(iter);
    fopt=inf;
    avgiter=0;
    yopt=[];
    for i=1:nbres
        if fval{i}<fopt
            fopt = fval{i};
            yopt=y{i};
        end
        avgiter = avgiter + iter{i};
    end
    avgiter=int32(avgiter/nbres);
end