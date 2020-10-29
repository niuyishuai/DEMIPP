function c = computec(grad,x,r)
H=jacobian(grad,x);
n=length(x);
L=zeros(1,n);
for i=1:n
    for j=1:n
        [cc,mm]=coefficients(H(i,j));
        nn=length(cc);
        for k=1:nn
            L(i)=L(i)+abs(cc(k))*r^degree(mm(k));
        end
    end
end
c=max(L);
end