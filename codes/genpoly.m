function [J,x,base,coef]=genpoly(N,d,polytype,density)
    % polytype: 0: polylab, 1: multipoly, 2: yalmip, 3: syms, 4: sostools
    
    % set default parameters
    if nargin < 4
        density = 1;
    end
    if nargin <3
        polytype = 1;
    end
    
    % create polynomial variables and basis
    switch polytype
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
    coef=sprand(nb,1,density);
    idx=[coef~=0];
    coef(idx)=round(20*coef(idx)-10);
    % for dense polynomial
    %coef=randi([-10,10],nb,1);
    
    % create polynomial
    switch polytype
        case 0
            J=simplify(coef'*base);
        case 1
            J = mpoly_mtimes_leftdoublematrix(coef',base);
            %J = mpoly_mtimes(coef',base);
            %J = mpoly_simplify(J);
        otherwise
            J=coef'*base; % polynomial objective function
    end
end