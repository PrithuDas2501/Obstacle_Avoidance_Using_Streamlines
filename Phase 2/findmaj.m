function y = findmaj(q,d,c,A)
a = ((q-c)'*A*(q-c)-1)^2;
b = norm(q-d);

y = a - b;
end