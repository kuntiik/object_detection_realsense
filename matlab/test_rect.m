

width = 10;
height = 30;
cx = 20;
cy = 40;
center = [cx; cy];

rect_whc = [cx - width/2, cy - height/2, width, height];
points = bbox2points(rect_whc)';
points = [points; ones(1,4)];


theta = pi/6;
rot = [cos(theta) -sin(theta) 0;
       sin(theta) cos(theta) 0;
        0 0 1];
trans_fw = [1 0 -cx;
            0 1 -cy;
            0 0  1];
trans_bw = [1 0 cx;
            0 1 cy;
            0 0 1];
rotate_around_center = trans_bw*rot*trans_fw;

points(:,end+1) = points(:,1);
points = rotate_around_center* points;

vec = points(:,1) - points(:,4);
norm = [vec(2); -vec(1)];
offs = norm'*center;
x_line = linspace(0,100);
syms a b c x y l x0 y0;
eq1 = a*x + b*y - c == 0;
eq2 = sqrt((x -x0)^2 + (y - y0)^2) == l;
sol = solve(eq1, y);
sol = subs(sol, [a b c], [norm(1) norm(2) offs]);
sol_y = double(subs(sol, x, x_line));

sol_pnt = solve([eq1,eq2],[x,y]);
sol_pntx = double(subs(sol_pnt.x, [a b c l x0 y0], [norm(1) norm(2) offs 20 cx cy]));
sol_pnty = double(subs(sol, x, sol_pntx));

hold on
plot(points(1,:),points(2,:), '*-');
plot(x_line, sol_y);
plot(sol_pntx, sol_pnty, 'o');
hold off


