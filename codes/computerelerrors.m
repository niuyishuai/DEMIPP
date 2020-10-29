count=0;
solutions=zeros(size(listofsols,1),4);
for N=setN
    for d=setd
        count=count+1;
        x = polylabvar(N);
        base = MPOLY.monolist(N,d);
        %x = sdpvar(N,1);
        %base=monolist(x,d);
        coef=listofprobs{count};
        J=coef'*base;
        
        % compute the values of J at best solutions
        solutions(count,1)=evalfcn(J,x,round(listofsols{count,1}));
        solutions(count,2)=evalfcn(J,x,round(listofsols{count,2}));
        solutions(count,3)=evalfcn(J,x,round(listofsols{count,3}));
        solutions(count,4)=evalfcn(J,x,round(listofsols{count,4}));        
    end
end
disp('Round solutions computed!');

%%
errors=nan(size(listf,1),4);
for i=1:size(listf,1)
    errors(i,1)=abs(solutions(i,1)-listf(i,5))/(1+abs(listf(i,5)));
    errors(i,2)=abs(solutions(i,2)-listf(i,5))/(1+abs(listf(i,5)));
    errors(i,3)=abs(solutions(i,3)-listf(i,5))/(1+abs(listf(i,5)));
    errors(i,4)=abs(solutions(i,4)-listf(i,5))/(1+abs(listf(i,5)));
end
disp('All errors computed!');
