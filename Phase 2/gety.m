function y = gety(Y,X,streamval,alpha,rx,ry)
A = Y*cos(alpha)/ry + X*sin(alpha)/rx;
B = 1/((X*cos(alpha)/rx - Y*sin(alpha)/ry)^2 + (Y*cos(alpha)/ry + X*sin(alpha)/rx)^2);

y = A*(1-B) - streamval;
end