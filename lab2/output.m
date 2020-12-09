x=infsup(-8,8);
y=infsup(-8,8);
z = globopt0([x,y]);
inters = [];

for i = 1:200
    xl = z(i).Box(1).inf;
    xr = z(i).Box(1).sup;
    yl = z(i).Box(2).inf;
    yr = z(i).Box(2).sup;
    rectangle('Position',[xl yl xr-xl yr-yl])
    %left = [z(i).Box(1).inf, z(i).Box(2).inf]
    %right = [z(i).Box(1).sup, z(i).Box(2).sup]
    inters = [inters; z(i).Box];
    % z(i).Box
end
inters
axis([-8 8 -8 8])
