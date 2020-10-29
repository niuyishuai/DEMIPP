function Feps = generatepenaltyfunc(J,x,c,epsilon)
    % compute penalized objective function
    if isa(J,'function_handle')
        Feps = @(xx)J(xx)+c*(xx'*xx)/2+(1/(4*epsilon))*(sum((xx.^2-1).^2));
    elseif isa(J,'struct')
        Feps = J;
        Feps = mpoly_plus(Feps, mpoly_mtimes_doublescalar(c/2, mpoly_mtimes_polymatrix(x',x)));
        w = mpoly_minus(mpoly_times(x,x),1);
        Feps = mpoly_plus(Feps, mpoly_mtimes_leftdoublematrix(1/(4*epsilon), mpoly_mtimes_polymatrix(w',w)));
        Feps = mpoly_simplify(Feps);
    else
        v = (x.^2-1);
        Feps = J+c*(x'*x)/2+(1/(4*epsilon))*(v'*v);
    end
end