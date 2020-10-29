function grad = computegradient(f,x)
    if isa(f,'sdpvar') || isa(f,'polynomial')
        grad = jacobian(f,x)';
        %grad = yalmip2mpoly(grad,x);
    elseif isa(f,'sym')
        grad = gradient(f,x);
    elseif isa(f,'MPOLY')
        grad = f.jacobian';
    elseif isa(f,'struct')
        grad = mpoly_jacobian(f)';
        %grad = mpoly_simplify(grad);     
    end
end