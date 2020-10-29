%% compare local solvers
hold on
drawtimelist(listtime(:,1),'-.s');
drawtimelist(listtime(:,2),'-.>');
drawtimelist(listtime(:,3),'-.o');
drawtimelist(listtime(:,4),'-.+');
legend('Houbolt','Lie','RK(4,5)','IPOPT','Location','northwest');
setupfig;
ylabel('Log accumulated computation time');
xlabel('Percentage of problems solved');
