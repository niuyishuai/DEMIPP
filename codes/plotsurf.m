function plotsurf(figid,F,x,xbounds,ybounds)
% figid: figure id (often take 1)
% F: function to be ploted
% x: function variable
% xbounds, ybounds: bounds of x and y, as [xmin, xmax], [ymin, ymax];
figure(figid);
view(3);
[X,Y]=meshgrid([xbounds(1):0.15:xbounds(2)],[ybounds(1):0.15:ybounds(2)]);
Z=zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j) = evalfcn(F,x,[X(i,j);Y(i,j)]);
    end
end
surf(X,Y,Z);
%xlabel('v_1');
%zlabel('\Pi');
%ylabel('v_2');
alpha(0.5); % set transparency
drawnow limitrate;
end