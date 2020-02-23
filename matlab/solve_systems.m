clear all;

syms x y a b c px py l; 
p = [-383, 412, 1612]';
brick_n = [0.464 0.886]';

offs = p(1:2)'*brick_n;

eq1 = a*x + b*y - c == 0;
eq2 = sqrt((x - xp)^2 + (y - yp)^2) == l;
% eq1 = brick_n(1)*x + brick_n(2)*y - offs == 0;
% eq2 = sqrt((x -p(1))^2 + (y - p(2))^2) == 300;
% x_axis = linspace(-440, -330);
% sol = solve(eq1, y);
% y_axis = double(subs(sol, x, x_axis));
% 
% sol2 = solve([eq1,eq2]);
% xp = double(sol2.x);
% yp = double(sol2.y);



hold on;
plot(x_axis,y_axis);
plot(p(1), p(2), 'bo', 'MarkerSize', 10);
plot(xp, yp, 'gx', 'MarkerSize', 6);
hold off


