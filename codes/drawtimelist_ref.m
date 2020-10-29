function [xx,yy]=drawtimelist(timelist,linetype)
    %figtype='avgeachn';
    figtype='logavgeachn';
    switch figtype
        case 'small2large'
            lst=small2large(timelist);
        case 'large2small'
            lst=large2small(timelist);
        case 'avgeachn'
            range_n=4:2:20;
            yy=mean(reshape(timelist,3,length(range_n)));
            hold all;
            plot(range_n,yy,linetype,'LineWidth',1.5);
            return;
        case 'logavgeachn'
            range_n=4:2:20;
            yy=log(mean(reshape(timelist,3,length(range_n))));
            hold all;
            plot(range_n,yy,linetype,'LineWidth',1.5);
            return;
    end
    
    n=length(timelist);
    xx=1:n;
    yy=zeros(n,1);
    yy(1)=lst(1);
    for i=2:n
        yy(i)=lst(i)+yy(i-1);
    end
    hold all
    plot(xx,yy,linetype,'LineWidth',1.5);
end

function lst=small2large(x)
    lst=sort(x);
end

function lst=large2small(x)
    lst=sort(x,'descend');
end