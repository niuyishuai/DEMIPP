function p=randonunitsphere(n)
    p=randn(n,1);
    p=p/norm(p);
end