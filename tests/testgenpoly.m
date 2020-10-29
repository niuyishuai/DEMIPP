Nspan=2:12;
dspan=2:4;
density=1; % density of polynomial
tlst=[];            
for N=Nspan
    for d=dspan
        for polytype=[0 2 3 4] % 1: multipoly, 2: yalmip, 3: syms, 4: sostool
            tic;
            [J,x]=genpoly(N,d,polytype,rand);
            x0=2*rand(N,1)-1;
            evalfcn(J,x,x0);
            t(polytype+1)=toc;
        end
        fprintf('%d x %d: %s\n',N,d,num2str(t));
        tlst=[tlst;t];
    end
end

%%
fprintf('mean time: %s\n',num2str(mean(tlst)));
fprintf('std time: %s\n',num2str(std(tlst)));
figure;
hold all; 
drawtimelist(tlst(:,1),'-o');
drawtimelist(tlst(:,3),'-s');
drawtimelist(tlst(:,4),'-+');
drawtimelist(tlst(:,5),'-.');
%legend('multipoly','sostool');
legend('polylab','yalmip','matlab sym','sostool'); %'matlab sym',
title(sprintf('N=[%d,%d], d=[%d,%d]',min(Nspan),max(Nspan),min(dspan),max(dspan)));
xlabel('number of problems');
ylabel('total time (sec.)');