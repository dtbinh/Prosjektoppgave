syms x x0 y y0 z z0

f = symfun(sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2), [x, x0, y, y0, z, z0]);

jaco = diff(f, x);