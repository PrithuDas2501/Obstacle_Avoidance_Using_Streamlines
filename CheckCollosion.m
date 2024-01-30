function y = CheckCollosion(lambda, A, c, Pos, Vel)

RelPos = Pos - c;

l = (RelPos + lambda*(Vel));

y = l'*A*l - 1;
end