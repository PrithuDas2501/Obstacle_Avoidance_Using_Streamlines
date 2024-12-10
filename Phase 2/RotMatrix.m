function rm = RotMatrix(Theta, Zdash)
eulzxy = quat2eul([cos(Theta/2) sin(Theta/2)*Zdash']);

% Thetax = Theta*Zdash(1);
% Thetay = Theta*Zdash(2);
% Thetaz = Theta*Zdash(3);

Thetax = eulzxy(3);
Thetay = eulzxy(2);
Thetaz = eulzxy(1);

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