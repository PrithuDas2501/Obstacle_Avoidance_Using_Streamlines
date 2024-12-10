function y = find_stream2(r,alpha,xdashdash,ydashdash,rx,ry)
X = dot(r,xdashdash);
Y = dot(r,ydashdash);

A = Y*cos(alpha)/ry + X*sin(alpha)/rx;
B = 1/((X*cos(alpha)/rx - Y*sin(alpha)/ry)^2 + (Y*cos(alpha)/ry + X*sin(alpha)/rx)^2);

y = A*(1-B);
y
end