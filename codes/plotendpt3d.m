function plotendpt3d(figid,x)
% plot ending point with red color
figure(figid);
scatter3(x(1),x(2),x(3),80,'filled','o','MarkerFaceColor','y');
drawnow limitrate
end