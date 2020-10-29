%% order of delta
%setN=[2 4 6 8 10 12 14 16];
%setd=[4 6];
%setN=[2 4 6 8 10];
%setd=[4 5 6];
%param.epsilon=1e-4; % 1e-3;
%param.c=10;
%param.m=1;
%param.gamma=50;

count=0;
figure;
hold on;
for N=setN
    for d=setd
        count=count+1;
        %phi(count)=param.epsilon*N*d;
        phi = sqrt(N)/param.epsilon;
        %scatter(log(listdelta(count,1)),phi,'ro');
        %scatter(log(listdelta(count,2)),phi,'bs');
        scatter(log(listdelta(count,3)),phi,'mo');
        xlabel('\delta');
        ylabel('Nd\epsilon');
    end
end