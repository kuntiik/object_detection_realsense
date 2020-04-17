clear all;
syms x y z x0 y0 z0 a b c d p m n k l; 

%symbolic equations
eq1 = a*x + b*y  + c*z - d == 0;
eq2 = p*x + m*y +n*z -k == 0;
eq3 = (x - x0)^2 + (y - y0)^2 +(z - z0)^2 == l^2;


%solve with respect to y
sol = solve([eq1,eq2, eq3],[x,y,z]);
x1 = simplify(sol.x(1));
x2 = simplify(sol.x(2));
y1 = simplify(sol.y(1));
y2 = simplify(sol.y(2));
z1 = simplify(sol.z(1));
z2 = simplify(sol.z(2));