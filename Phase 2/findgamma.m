function g = findgamma(A,c,x)
x = x-c;
g = x'*A*x-1;
end