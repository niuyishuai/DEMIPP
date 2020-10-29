function [xx,yy]=drawtimelist(timelist,linetype)
n=length(timelist);
xx=[0:n-1]*100/(n-1);
yy=zeros(1,n);
yy(1)=timelist(1);
for i=2:n
    yy(i)=timelist(i)+yy(i-1);
end
xlim([0,100]);
plot(xx,log(yy),linetype,'LineWidth',1.5);
end