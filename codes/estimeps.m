function epsilon = estimeps(zeta,c,J,x)
    % given the precision zeta, estimate an upper bound for epsilon
    
    N=getnbvars(J,x);
    grad=computegradient(J,x);
    for i=1:length(grad)
        grad(i).coef=abs(grad(i).coef);
    end
    r=sqrt(N);
    Lr=max(grad.eval(r*ones(N,1)));
    epsilon = zeta/(4*(c+Lr)*r);    
end