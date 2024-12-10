clear all; close all; clc
%%
a = 1;
r = -10*a:0.1:10*a;
rx = 5;
ry = 2;
[x,y] = meshgrid(r,r);
d = -pi/4;
%theta = linspace(0,2*pi,100);
e = cos(d) + 1i*sin(d);
z = (x/rx + 1i*y/ry)*e;
%R = real(z)/rx;
%I = imag(z)/ry;
%z = R+i*I;
psi = imag(z+1./z);
%M=contour(r,r,psi,[0, 0.25, 0.5, 0.75, 1, -0.25, -0.5, -0.75, -1],'r');%% Keep as [1,-1] for best path
M=contourc(r,r,psi,[0.3,0.3]);
plot(M(1,2:1+M(2,1)), M(2,2:1+M(2,1)),'r');

%%
[x,y,z] = C2xyz(M);
%%
figure
hold on
grid on
axis equal
plot3(l(1,2,:)', l(2,2,:)', l(3,2,:)','r');
%plot(M(1,2:1+M(2,1)), M(2,2:1+M(2,1)),'r');
%plot(M(1,end/2+2:end/2 + M(2,end/2+1)), M(2,end/2+2:end/2 + M(2,end/2+1)),'r');
%plot(rx*a*cos(theta),ry*a*sin(theta))
%plot3(M(1,:),M(2,:),M(3,:),'b')
%%
x0 = 14;
y0 = 8 ;
z0 = 0 ;
A = [1/4 0 0; 0 1/9 0; 0 0 1/25];
%% [Aye,Bye,qx,qy,qz]=intersection_elipsoid_plane2([A(1,1) A(2,2) A(3,3) (A(2,1)+A(1,2))/2 (A(3,1)+A(1,3))/2 (A(2,3)+A(3,2))/2 (-2*A(1,1)*x0-(A(2,1)+A(1,2))*y0-(A(3,1)+A(1,3))*z0)/2 (-2*A(2,2)*y0-(A(2,1)+A(1,2))*x0-(A(3,2)+A(2,3))*z0)/2 (-2*A(3,3)*z0-(A(3,2)+A(2,3))*y0-(A(3,1)+A(1,3))*x0)/2 A(1,1)*x0^2+A(2,2)*y0^2+A(3,3)*z0^2 + (A(2,1)+A(1,2))*x0*y0 + (A(3,1)+A(1,3))*z0*x0 + (A(2,3)+A(3,2))*y0*z0-1], [1 1 1 1])
ra = svd(A);
Zdash = [0;0;1];
plane = [Zdash' -4];
for q = 1:4
    if abs(plane(q))<0.001
        plane(q) = 0.001;
    end
end
%[Aye,Bye,qx,qy,qz]=intersection_elipsoid_plane2([A(1,1) A(2,2) A(3,3) (A(2,1)+A(1,2))/2 (A(3,1)+A(1,3))/2 (A(2,3)+A(3,2))/2 (-2*A(1,1)*x0-(A(2,1)+A(1,2))*y0-(A(3,1)+A(1,3))*z0)/2 (-2*A(2,2)*y0-(A(2,1)+A(1,2))*x0-(A(3,2)+A(2,3))*z0)/2 (-2*A(3,3)*z0-(A(3,2)+A(2,3))*y0-(A(3,1)+A(1,3))*x0)/2 A(1,1)*x0^2+A(2,2)*y0^2+A(3,3)*z0^2 + (A(2,1)+A(1,2))*x0*y0 + (A(3,1)+A(1,3))*z0*x0 + (A(2,3)+A(3,2))*y0*z0-1], [1 1 1 1])

[Aye,Bye,qx,qy,qz]=intersection_elipsoid_plane2([A(1,1) A(2,2) A(3,3) (A(2,1)+A(1,2))/2 (A(3,1)+A(1,3))/2 (A(2,3)+A(3,2))/2 0 0 0 -1], plane)
Oe = [0;0;0];
Xg = [5;5;4];
Xdash = [qx;qy;qz] - Xg + Oe;
Xdash = Xdash/norm(Xdash);
Zdash = Zdash/norm(Zdash);
Ydash = cross(Zdash,Xdash);
%%%%%%% Add the fmincon function for finding the angle of attack
%%
AOA = fmincon(@(theta) FindAOA(theta, Aye, [qx;qy;qz], Xdash, Zdash, A),0.1, [1;-1], [pi/2;pi/2]);
Non_Transformed_Plane_AOA = atan((Aye/Bye)*tan(AOA));
%%
a = 1;
r = -10*a:0.1:10*a;
[x,y] = meshgrid(r,r);
theta = linspace(0,2*pi,100);
e = cos(Non_Transformed_Plane_AOA) + 1i*sin(Non_Transformed_Plane_AOA);
z = (x/Aye + 1i*y/Bye)*e;
%R = real(z)/rx;
%I = imag(z)/ry;
%z = R+i*I;
psi = imag(z+1./z);
figure
hold on
grid on
axis equal
M=contourc(r,r,psi,[0.25,-0.25]); %% Keep as [1,-1] for best path
plot(M(1,2:1+M(2,1)), M(2,2:1+M(2,1)),'r');
plot(M(1,end/2+2:end/2 + M(2,end/2+1)), M(2,end/2+2:end/2 + M(2,end/2+1)),'r');
plot(Aye*a*cos(theta),ry*Bye*sin(theta))
%%
Rotation_Matrix = RotMatrix(AOA, Zdash);
Xdashdash = Rotation_Matrix*Xdash;
Ydashdash = Rotation_Matrix*Ydash;
figure
hold on
grid on
ellipsoid(0,0,0,ra(1)^-0.5,ra(2)^-0.5,ra(3)^-0.5)
for i = 1:M(2,1)
    P = Xdashdash*M(1,1+i) + Ydashdash*M(2,1+i);
    plot3(P(1), P(2), P(3),'.-');
end
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
%{
% Define the plane equation: z = ax + by + c
a = -plane(1)/plane(3);
b = -plane(2)/plane(3);
c = -plane(4)/plane(3);

% Define the range of x and y values
[x, y] = meshgrid(-10:0.1:10, -10:0.1:10);

% Calculate the corresponding z values
z = a * x + b * y + c;

% Plot the plane using surf
surf(x, y, z, 'FaceAlpha', 0.7);
title('Plane in 3D');
%}
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

%%
function y = FindAOA(theta, Rx, C, Xdash, Zdash, A)
Rotation_Matrix = RotMatrix(theta,Zdash);
y = ((Rx*Rotation_Matrix*Xdash)'*A*(Rx*Rotation_Matrix*Xdash)-1)^2;
end

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
%{
theta = Got from the fmincon
Find the stream function value
Follow stream line for this value of stream function
%}