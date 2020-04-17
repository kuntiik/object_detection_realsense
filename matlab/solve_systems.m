clear all;
%i have center of brick represented by point in 3d and vector of plane and
%brick direction, need to compute rectangle edges
syms x y z a b c l x0 y0 z0 k m n p; 
%some sample point and vetor

%symbolic equations
eq1 = a*x + b*y  - c == 0;
%eq2 = sqrt((x - x0)^2 + (y - y0)^2) == l;
eq2 = (x - x0)^2 + (y - y0)^2 == l^2;
%eq2_2 = sqrt((x - x0)^2 + (y - y0)^2 +(z - z0)^2) == l;
eq2_2 = (x - x0)^2 + (y - y0)^2 +(z - z0)^2 == l^2;
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

%corners = [119.336, 13.1241, 1012.28;
          %-30.9463, 109.701, 1102.21;
          %-301.316, -311.017, 1433.72;
          %-151.034, -407.594, 1343.79]';
%pl_n = [0.600948, 0.219414, 0.768583]';
%br_n = [-0.841262, 0.540627]';

corners = [-675.018, 204.22, 1640.49;
-518.223, 209.321, 1516.44;
-537.102, 789.719, 1365.51;
-693.898, 784.618, 1489.56]';
pl_n = [0.600948, 0.219414, 0.768583]';
br_n = [0.999471, 0.0325106]';
cor_n = corners(:,1) - corners(:,2);
cor_n = cor_n/norm(cor_n);
cor_n_2d = cor_n(1:2,:)/norm(cor_n(1:2,:));

c_pl = rad2deg(acos(dot(cor_n, pl_n)))
c_br = rad2deg(acos(dot(cor_n_2d, br_n)))

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


