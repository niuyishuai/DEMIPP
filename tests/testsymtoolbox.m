% create polynomial using symbolic toolbox
N=2;
d=4;
x=sym('x',[N,1],'real');
%syms 'x%d' [N 1] real
basis=monolist(x,d);
nb=length(basis);
coef=randi([-10,10],nb,1);
J=coef'*basis;

% check degree of polynomial
polynomialDegree(J)

% get coefficients and monomials
[cc,bb]=coeffs(J)

% eval polynomial
subs(J,x,[1;2]) % Return a symbolic result
double(subs(J,x,[1;2]))

% compute gradient, jacobian and hessian
gradient(J,x)
jacobian(J,x)
H=hessian(J,x)

% compute norms
norm(H,1)
norm(H,2)
norm(H,inf)
norm(H,'fro')
