function nbvars = getnbvars(f,x)
    if isa(f,'sdpvar') || isa(f,'polynomial') || isa(f,'sym')
        nbvars = numel(x);
    elseif isa(f,'MPOLY') || isa(f,'struct')
        nbvars = f.n;
    end
end