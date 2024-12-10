function rm = RotMatrix(Theta, Zdash)
Thetax = Theta*Zdash(1);
Thetay = Theta*Zdash(2);
Thetaz = Theta*Zdash(3);

Rx = [1 0           0;
      0 cos(Thetax) -sin(Thetax);
      0 sin(Thetax) cos(Thetax)];

Ry = [cos(Thetay)  0 sin(Thetay);
      0            1 0;
      -sin(Thetay) 0 cos(Thetay)];

Rz = [cos(Thetaz) -sin(Thetaz) 0;
      sin(Thetaz) cos(Thetaz)  0;
      0           0            1];

rm = Rx*Ry*Rz;
end