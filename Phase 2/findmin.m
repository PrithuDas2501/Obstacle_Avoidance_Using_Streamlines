function y = findmin(q,d,c,A)

options = optimset('Display','off');
k = fsolve(@(k) ((d + k*(q-d) - c)'*A*(d + k*(q-d) - c)-1)^2,1,options);
b = (norm(k*(q-d)))^2;

y = b;
end