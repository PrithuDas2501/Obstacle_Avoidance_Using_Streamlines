function y = elli(x,xc,yc,a,b)

f = @(x,ycord) ((x-xc)/a)^2 + ((ycord-yc)/b)^2 - 1;

options = optimset('Display','off');
y = fsolve(@(z) f(x,z), b+yc, options);

end