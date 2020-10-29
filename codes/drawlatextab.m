%% format 1
% N | d | Houbolt | Lie | RK(4,5)
idx=1;
for N=setN
    for d=setd
        fprintf('$%d$ & $%d$ & ',N,d);
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,1),listiter(idx,1),listtime(idx,1),listdelta(idx,1));  
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,2),listiter(idx,2),listtime(idx,2),listdelta(idx,2));
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ \\\\\n',listf(idx,3),listiter(idx,3),listtime(idx,3),listdelta(idx,3));
        idx=idx+1;
    end
end
fprintf('\\hline\n');
fprintf('\\multicolumn{2}{c|}{average} & & $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,1))),mean(listtime(:,1)),mean(listdelta(:,1)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,2))),mean(listtime(:,2)),mean(listdelta(:,2)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ \n',int32(mean(listiter(:,3))),mean(listtime(:,3)),mean(listdelta(:,3)));

%% format 2
% N | d | Houbolt | Lie | RK(4,5) | Exhaustive
idx=1;
for N=setN
    for d=setd
        fprintf('$%d$ & $%d$ & ',N,d);
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,1),listiter(idx,1),listtime(idx,1),listdelta(idx,1));  
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,2),listiter(idx,2),listtime(idx,2),listdelta(idx,2));
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,3),listiter(idx,3),listtime(idx,3),listdelta(idx,3));
        fprintf('$%.3e$ & $%.2f$ \\\\\n',listf(idx,4),listtime(idx,4));
        idx=idx+1;
    end
end
fprintf('\\hline\n');
fprintf('\\multicolumn{2}{c|}{average} & & $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,1))),mean(listtime(:,1)),mean(listdelta(:,1)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,2))),mean(listtime(:,2)),mean(listdelta(:,2)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ \n',int32(mean(listiter(:,3))),mean(listtime(:,3)),mean(listdelta(:,3)));
fprintf('& $%.2f$ \n',mean(listtime(:,4)));

%% format 3
% N | d | Houbolt | Lie | RK(4,5) | QuadReduction + Gurobi | Exhaustive
idx=1;
for N=setN
    for d=setd
        fprintf('$%d$ & $%d$ & ',N,d);
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,1),listiter(idx,1),listtime(idx,1),listdelta(idx,1));  
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,2),listiter(idx,2),listtime(idx,2),listdelta(idx,2));
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,3),listiter(idx,3),listtime(idx,3),listdelta(idx,3));
        fprintf('$%.3e$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,4),listiter(idx,4),listtime(idx,4),listdelta(idx,4));  
        fprintf('$%.3e$ & $%.2f$ \\\\\n',listf(idx,5),listtime(idx,5));
        idx=idx+1;
    end
end
fprintf('\\hline\n');
fprintf('\\multicolumn{2}{c|}{average} & & $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,1))),mean(listtime(:,1)),mean(listdelta(:,1)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,2))),mean(listtime(:,2)),mean(listdelta(:,2)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ \n',int32(mean(listiter(:,3))),mean(listtime(:,3)),mean(listdelta(:,3)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ \n',int32(mean(listiter(:,4))),mean(listtime(:,4)),mean(listdelta(:,4)));
fprintf('& $%.2f$ \n',mean(listtime(:,4)));

%% format 5
% N | d | Houbolt | Lie | RK(4,5) | Ipopt | QuadReduction + Gurobi
idx=1;
for N=setN
    for d=setd
        fprintf('$%d$ & $%d$ & ',N,d);
        fprintf('$%d$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,1),listiter(idx,1),listtime(idx,1),listdelta(idx,1));  
        fprintf('$%d$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,2),listiter(idx,2),listtime(idx,2),listdelta(idx,2));
        fprintf('$%d$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,3),listiter(idx,3),listtime(idx,3),listdelta(idx,3));
        fprintf('$%d$ & $%d$ & $%.2f$ & $%.2e$ &',listf(idx,4),listiter(idx,4),listtime(idx,4),listdelta(idx,4));
        fprintf('$%d$ & $%.2f$ \\\\\n',listf(idx,5),listtime(idx,5)); 
        idx=idx+1;
    end
end
fprintf('\\hline\n');
fprintf('\\multicolumn{2}{c|}{average} & & $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,1))),mean(listtime(:,1)),mean(listdelta(:,1)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,2))),mean(listtime(:,2)),mean(listdelta(:,2)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ & ',int32(mean(listiter(:,3))),mean(listtime(:,3)),mean(listdelta(:,3)));
fprintf('& $%d$ & $%.2f$ & $%.2e$ ',int32(mean(listiter(:,4))),mean(listtime(:,4)),mean(listdelta(:,4)));
fprintf('& & $%.2f$ \\\\\n',mean(listtime(:,5)));

%% draw table of relative errors
% N | d | Houbolt | Lie | RK(4,5) | QuadReduction + Gurobi
idx=1;
for N=setN
    for d=setd
        fprintf('$%d$ & $%d$ & ',N,d);
        fprintf('$%.2f$ &',errors(idx,1));
        fprintf('$%.2f$ &',errors(idx,2));
        fprintf('$%.2f$ &',errors(idx,3));
        fprintf('$%.2f$ \\\\\n',errors(idx,4));
        idx=idx+1;
    end
end
fprintf('\\hline\n');
mm=mean(errors);
fprintf('\\multicolumn{2}{c|}{average} & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \n',mm(1),mm(2),mm(3),mm(4));
