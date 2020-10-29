function [xopt,fval,iter,status] = callgurobimip(Q,q,c,y0)
% Call Gurobi to solve boolean quadratic optimization
%  minimize
%        x'*Q*x + q'*x + c
%  subject to
%        x binaries
xopt=[];
fval=[];
status=[];
N = length(q);
model.vtype = 'B';
model.modelsense = 'min';

names=cell(N,1);
for i=1:N
    names{i}=sprintf('x%d',i);
end
model.varnames = names;

model.Q = Q;
model.obj = full(q);
model.start = y0;
model.A = sparse(1,N);
model.b = 0;

%gurobi_write(model, 'bqp.lp');

params.outputflag = 0;

result = gurobi(model, params);
%disp(result);

if strcmp(result.status,'OPTIMAL')
    xopt=result.x;
    fval=result.objval+c;
    iter=result.itercount;
    status=result.status;
end
end
