function y = Inplane(X, Xdashdash, Ydashdash, proj)
y = X(1)*Xdashdash + X(2)*Ydashdash - proj;
end