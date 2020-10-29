function [f,Q,q,c]=reducemin(P,x,polytype)
    %reducemin
    %% DESCRIPTION:
    %  Reduce an n variate polynomial P(x) over {0,1}^n as a quadratic
    %  polynomial z'*Q*z+q'z+c.
    %
    %% SYNTAX:
    %  p = reducemin(p,[],0);
    %
    %% INPUTS:
    %  P: polynomial of polytype
    %  x: polynomial variables (unused)
    %  polytype: 0: polylab (class) 1: polylab (function), 2: yalmip, 3: syms, 4: sostools
    %
    %% OUTPUTS:
    %  Q,q,c: coefficients of the quadratic polynomial z'*Q*z+q'z+c
    %
    %% EXAMPLE:
    %   p = MPOLY(3,[1;2;3],[2 0 0; 1 2 3; 0 1 1]);
    %   [Q,q,c]=reducemin(p,[],0)
    %
    %%
    %  See also MPOLY, mpvar, sym, sdpvar
    %
    %% COPYRIGHT:
    %  Copyright 2020, Yi-Shuai NIU. All Rights Reserved.
    %  2020/08/21    Initial Coding
    
    % extract data of polynomial
    switch polytype
        case 0 % polylab (class)
            n=P.n;
            coef=P.coef;
            pow=P.pow;
        case 1 % polylab (function)
        case 2 % yalmip
        case 3 % Matlab sym
        case 4 % sostools
            n=P.nvars;
            coef=P.coefficient;
            pow=P.degmat;
    end
    
    % reduce degree x^k terms to x (since x^k = x for x in {0,1})
    pow=double(pow>0);
    
    % combine common terms (it is possible since x^k is replaced by x)
    [pow,coef]=combinecommonterms(pow,coef);
    
    % reduce to quadratic
    M =1 + sum(abs(coef)); % 1 + 2*sum(abs(coef));
    m=n;
    
    %[pow,coef]=simplifypoly(pow,coef);
    [maxdeg,idx]=max(sum(pow,2));
    while maxdeg>2
        row=pow(idx,:);
        S=find(row);
        i=S(1);
        j=S(2);
        [eqpos,includepos]=locatesubsetpos(pow,[i,j]);
        if isempty(coef(eqpos)) % it is possible that no [i,j] term alone
            coef(end+1)=M;
            pow(end+1,m+1)=1;
            pow(end,:)=gennewrow(m+1,[i,j]);
        else % modify [i,j] term
            coef(eqpos)=coef(eqpos)+M;
        end
        coef(end+1)=-2*M;
        coef(end+1)=-2*M;
        coef(end+1)=3*M;
        pow(end+1,m+1)=1;
        pow(end,:)=gennewrow(m+1,[i,m+1]);
        pow(end+1,:)=gennewrow(m+1,[j,m+1]);
        pow(end+1,:)=gennewrow(m+1,[m+1]);
        nbadded=length(includepos);
        for ii=1:nbadded
            coef(end+1)=coef(includepos(ii));
            pow(end+1,m+1)=1;
            pow(end,:)=gennewrow(m+1,[setdiff(S,[i,j]),m+1]);
        end
        coef(includepos)=[];
        pow(includepos,:)=[];
        m=m+1;
        %[pow,coef]=simplifypoly(pow,coef);
        [maxdeg,idx]=max(sum(pow,2));
    end
    f=generatepoly(maxdeg,coef,pow,polytype);
    [Q,q,c]=formquad(pow,coef);
end

function [eqpos,includepos]=locatesubsetpos(pow,S)
    % find locations of all terms containing S, where eqpos for index of terms equals to S
    % and includepos for index of terms includes S.
    n=size(pow,2);
    v=zeros(1,n);
    v(S)=1;
    D=pow-v;
    % eqpos (position equals to S)
    eqpos=find(sum(abs(D),2)==0);
    % get includepos (position includes S)
    includepos=find([sum(D<0,2) + prod(D==0,2)]==0);
end

function row=gennewrow(len,S)
    % generate a new row with 1 for elements of indexes in S and 0 for others
    row=zeros(1,len);
    row(S)=1;
end

function [pow,coef]=simplfypoly(pow,coef)
    % eliminate all terms with zero coefficients
    idx=abs(coef)<1e-8;
    coef(idx)=[];
    pow(idx,:)=[];
end

function [pow,coef]=combinecommonterms(tpow,tcoef)
    % combine all common terms of polynomial
    
    % Sort pow
    [tpow,sortidx]= sortrows(tpow);
    % in the sorted power, find non repeated rows
    repidx=[1;any(tpow(2:end,:)~=tpow(1:end-1,:),2)];
    % extract non repeated pow
    tpow=tpow(repidx==1,:);
    % sum repeated coefs (using sparse matrix)
    tcoef=sparse(cumsum(repidx),sortidx,1)*tcoef;
    
    % Eliminate zero coefs
    idx = find(abs(tcoef)>1e-8);
    pow=tpow(idx,:);
    coef=tcoef(idx);
end

function [Q,q,c]=formquad(pow,coef)
    % formulate quadratic function with pow and coef
    n=size(pow,2);
    Q=zeros(n,n);
    q=zeros(n,1);
    % get constant term
    SM=sum(pow,2);
    idx=[SM==0];
    c=coef(idx);
    if isempty(c)
        c=0;
    end
    % get linear terms
    idx1=[SM==1];
    pow1=pow(idx1,:);
    coef1=coef(idx1);
    q(pow1*[1:n]')=coef1;
    % get quadratic terms
    idx2=[SM==2];
    pow2=pow(idx2,:);
    coef2=coef(idx2);
    for i=1:length(coef2)
        r=pow2(i,:);
        varpos=find(r>0);
        if length(varpos)==2
            Q(varpos(1),varpos(2))=coef2(i)/2;
            Q(varpos(2),varpos(1))=Q(varpos(1),varpos(2));
        else
            Q(varpos,varpos)=coef2(i);
        end
    end
    Q=sparse(Q);
    q=sparse(q);
end

function f=generatepoly(maxdeg,coef,pow,polytype)
    % generate a polynomial object from data
    switch polytype
        case 0 % polylab (class)
            f=MPOLY(maxdeg,coef,pow);
        case 1 % polylab (function)
        case 2 % yalmip
        case 3 % Matlab sym
        case 4 % sostools
    end
    
end
