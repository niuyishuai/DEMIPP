function Lr = estimLipparam(J,x,r)
    grad = computegradient(J,x);
    n=J.n;
    for i=1:n
        grad(i).coef=abs(grad(i).coef);
    end
    Lr = max(grad.eval(r*ones(n,1)));
end