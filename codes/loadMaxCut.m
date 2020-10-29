function  [J,x,Q,q,s,W]=loadMaxCut(filename)
    % read MaxCut problem from MQLip datafile. 
    % the problem is defined as
    % min 1/2 * sum_{i<j} W(i,j)*x(i)*x(j) - 1/2* sum_{i<j}W(i,j)
    % we return the matrix format (Q,q,s) for the quadratic objective function as:
    % 1/2*x'*Q*x + q'*x + s
    % The polylab format J,x and the coefficient matrix W is returned too.
    fid=fopen(filename,'r');
    getsize=false;
    currow=1;
    J=[];
    x=[];
    while ~feof(fid)
        tline = fgetl(fid);
        if tline(1)=='#' % skip comments
            continue
        end
        if getsize==false
            dim=sscanf(tline,'%d',[1,2]); % read problem size
            n=dim(1); % n number of nodes 
            m=dim(2); % m number of edges
            W=zeros(n,n); % adjacency matrix W
            %x=MPOLY.mpolyvars(n);
            k=m+1; % add 1 for constant term
            %pow=sparse(k,n);
            coef=zeros(k,1);
            getsize=true;
        else
            idx=sscanf(tline,'%f',[1,3]);
            i=idx(1);
            j=idx(2);
            v=idx(3);
            %pow(currow,i)=1;
            %pow(currow,j)=1;
            coef(currow)=v/2;
            W(i,j)=v;
            currow=currow+1;
        end
    end
    coef(end)=-sum(coef);
    %J=MPOLY(n,coef,pow);
    % symmetrize W
    W = W + W'; % the diagonal ara all zeros for MaxCut
    % generate matrix format
    Q=sparse(W/2);
    q=sparse(n,1);
    s=coef(end);
    fclose(fid);
end