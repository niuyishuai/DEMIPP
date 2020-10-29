function plotendpt2d(figid,x)
% plot ending point with red color
figure(figid);
scatter(x(1),x(2),80,'filled','o','MarkerFaceColor','y');
drawnow limitrate
end