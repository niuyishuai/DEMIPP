function vv = evalfcn(f,x,xx)
% evaluate the value of function of x at given point xx
if isa(f,'sdpvar')
    vv=replace(f,x,xx);
elseif isa(f,'sym')
    vv=double(subs(f,x,xx));
elseif isa(f,'struct')
    vv=mpoly_eval(f,xx);
elseif isa(f,'MPOLY')
    vv=f.eval(xx);
elseif isa(f,'polynomial')
    vv=double(subs(f,x,xx));
elseif isa(f,'function_handle')
    vv=f(xx);
else
    error('no supported function/variable type');
end
end

function fval = mysubs(f,x,xx)
fval=zeros(size(f));
for idx=1:numel(f)
    [coef,monos]=coeffs(f(idx),x);
    k=length(coef); % number of monomials
    n=length(x); % number of variables
    pow=zeros(n,1);
    for i=1:k % the i-th monomial
        for j=1:n % variable j
            pow(j)=polynomialDegree(monos(i),x(j));
        end
        fval(idx)=fval(idx)+coef(i)*prod(xx(:).^pow);
    end
end
end