function y = find_stream(r,alpha,ALPHA,rx,ry)
X = -r*cos(ALPHA);
Y = r*sin(ALPHA);

A = Y*cos(alpha)/ry + X*sin(alpha)/rx;
B = 1/((X*cos(alpha)/rx - Y*sin(alpha)/ry)^2 + (Y*cos(alpha)/ry + X*sin(alpha)/rx)^2);

y = A*(1-B);
end