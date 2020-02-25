  clear all;
%i have center of brick represented by point in 3d and vector of plane and
%brick direction, need to compute rectangle edges
syms x y z a b c l x0 y0 z0 k m n p; 
%some sample point and vetor

%symbolic equations
eq1 = a*x + b*y + 0*z - c == 0;
eq2 = sqrt((x - x0)^2 + (y - y0)^2) == l;
eq2_2 = sqrt((x - x0)^2 + (y - y0)^2 +(z - z0)^2) == l;
eq3 = p*x + m*y +n*z -k == 0;

% eq1 = brick_n(1)*x + brick_n(2)*y - offs == 0;
% eq2 = sqrt((x -p(1))^2 + (y - p(2))^2) == 300;
%x_axis = linspace(-440, -330);
%y_axis = double(subs(sol, x, x_axis));
% 

%solve with respect to y
sol = solve([eq1,eq2],[x,y]);
sol_3d = solve([eq1,eq2_2, eq3],[x,y,z]);

pt = [-383, 412, 1612]';
brick_n = [0.464 0.886]';
offs = pt(1:2)'*brick_n;


yp1 = double(subs(sol.y(1), [a,b,c,x0,y0,l],[brick_n(1), brick_n(2), offs, pt(1),pt(2),300]));
yp2 = double(subs(sol.y(2), [a,b,c,x0,y0,l],[brick_n(1), brick_n(2), offs, pt(1),pt(2),300]));

%i have y, now get x from equation of brick direction
sol_x = solve(eq1,x);
sol_x_num1 = double(subs(sol_x, [a,b,c,y],[brick_n(1), brick_n(2), offs, yp1]));
sol_x_num2 = double(subs(sol_x, [a,b,c,y],[brick_n(1), brick_n(2), offs, yp2]));

% 
% 
% hold on;
% plot(x_axis,y_axis);
% plot(p(1), p(2), 'bo', 'MarkerSize', 10);
% plot(xp, yp, 'gx', 'MarkerSize', 6);
% hold off


