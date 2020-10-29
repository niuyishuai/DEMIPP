function stamp = createtimestamp()
    l=fix(clock);
    str=[];
    for i=1:length(l-1)
        str=[str, num2str(l(i)),'_'];
    end
    str=[str,num2str(l(end))];
    stamp=str;
end