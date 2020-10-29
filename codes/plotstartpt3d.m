function plotstartpt3d(figid,x)
% plot starting point with yellow color
figure(figid);
scatter3(x(1),x(2),x(3),80,'filled','o','MarkerFaceColor','b');
drawnow limitrate
end