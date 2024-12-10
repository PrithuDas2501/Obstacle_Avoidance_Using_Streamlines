

% Arbitrary ELÝPSOÝD and A PLANE INTERSECTÝON


function[Aye,Bye,qx,qy,qz]=intersection_elipsoid_plane2(elp,plane)
%%%%%%
 
% Plane equation         AA.x + BB.y + CC.z+ DD = 0
% Standart Ellipsoid equation    (x/a)^2 + (y/b)^2 + (z/c)^2 = 1  
% Arbitrary Ellipsoid equation  Ax^2 + By^2 + Cz^2 +2 Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0


% Author: Sebahattin Bektas,Ondokuz Mayis University,2015 
%
% Inputs -----------------------------------------------------------------
%   Plane=[AA  BB  CC DD]   :Plane equation's coefficients  1x4
%                
% for standard ellipsoid    elp must be  1x3  =[a b c]
% for arbitrary ellipsoid   elp must be  1x10 =[ A  B C  D E  F G  H  I  J]
%
%   elp=[a b c]: Semi-axes of standart ellipsoid 
%    or
%   elp=[ A  B C  D E  F G  H  I  J]: Arbitrary Ellipsoid coefficient
% Outputs ----------------------------------------------------------------
%   Aye:  Semi-major axis of intersection ellipse
%   Bye:  Semi-minor axis of intersection ellipse
%
% qx,qy,qz : Cartesian coordinates of intersection ellipse's center
% This source code is given for free! However, I would be grateful if you refer 
% to corresponding article in any academic publication that uses 
% this code or part of it. 

% Please refer to: 

%  BEKTAS, S Orthogonal distance from an ellipsoid. Bol. Ciênc. Geod. [online]. 2014, vol.20, n.4, pp. 970-983. ISSN 1982-2170.

%  BEKTAS, S Intersection of an Ellipsoid and a Plane,International Journal of Research in Engineering and Applied Sciences VOLUME 6,ISSUE 6,2016 


%clc,clear all
  
 format long


   fprintf(' INTERSECTION ELÝPSOÝD and A PLANE \n by Sebahattin Bektas\n ===========================================\n  \n\n')
plane(1)=plane(1)/plane(4);plane(2)=plane(2)/plane(4);plane(3)=plane(3)/plane(4);plane(4)=plane(4)/plane(4);

Ax=plane(1);Ay=plane(2);Az=plane(3);

if length(elp)>3
   disp(' Arbitrary Ellipsoid equation  Ax^2 + By^2 + Cz^2 +2 Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0')
   fprintf('Arbitrary Ellipsoid equation A,B,C,D,E,F,G,H,I,J %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n',elp)
      fprintf('Plane equation  AA.x + BB.y + CC.z+ DD =0 \n')
      fprintf('Plane equation  AA BB CC DD %10.7f %10.7f %10.7f %10.7f\n',plane)

   %elp,plane
 %[center,radii,R,Angle,code]=Conic_EllipsoidParameter(elp)
 [center,radii,quat,R] = ellipsoid_im2ex(elp);
R=R'; R(:,2)=-R(:,2);
Tmat=[R center;0 0 0 1];

%düz
xyz=[1 2;2 3;3 1];

for i=1:3
    xyz(i,3)=(-1-Ax*xyz(i,1)-Ay*xyz(i,2))/Az;
end
xyz;
[param,artik] = duzdenkbelirleme(xyz);
 D21=[xyz ones(3,1)];
 kont=(inv(Tmat)*D21')';
 kont(:,4)=[];
    [param,artik] = duzdenkbelirleme(kont);

a=radii(1);b=radii(2);c=radii(3);
end
%%%%%%%%%%%%%%%% param=plane;
if length(elp)==3 
    param=plane;
 a=elp(1);b=elp(2);c=elp(3);   
  disp(' Standart Ellipsoid equation  (x/a)^2 + (y/b)^2 + (z/c)^2 = 1  ')
   fprintf('Standart Ellipsoid equation   a    b   c %10.7f %10.7f %10.7f \n',elp)
      fprintf('Plane equation  AA.x + BB.y + CC.z+ DD =0 \n')
      fprintf('Plane equation  AA BB CC  DD %10.7f %10.7f %10.7f %10.7f \n',plane)
end

%[Aye,Bye,q1,q2,q3]=EllipsoidPlaneIntersection(param(1),param(2),param(3),1,a,b,c);
%param;
A=param(1);B=param(2);C=param(3);D=1;
car=A*B*C*D;
if car~=0
kx2=1/a^2 + (A^2)/(C^2*c^2);
ky2=1/b^2 + (B^2)/(C^2*c^2);
kxy=(2*A*B)/(C^2*c^2);
kx=+ (2*A*D)/(C^2*c^2);
ky=+ (2*B*D)/(C^2*c^2) ;
ksab=D^2/(C^2*c^2)- 1;
ParA=[kx2  kxy  ky2  kx  ky  ksab];
G=AtoG(ParA);

q1=G(1);q2=G(2);q3=(A*q1+B*q2+D)/-C;
%q1,q2,q3
end

if C==0 & A~=0
ParA=[(1/b^2 + (B/A/a)^2) 0  1/c^2 (2*D*B/(A*a)^2)  0  (-1+(D/A/a)^2)];G=AtoG(ParA);

q2=G(1);q3=G(2);q1=(D+B*q2+C*q3)/-A   ;
end

if C==0 & B~=0
ParA=[(1/a^2 + (A/B/b)^2) 0  1/c^2 (2*D*A/(B*b)^2)  0  (-1+(D/B/b)^2)];G=AtoG(ParA);

q1=G(1);q3=G(2);q2=(D+A*q1+C*q3)/-B;    
end

if A==0 & B==0 ,    q1=0;q2=0;q3=-D/C;end 
if D==0,    q1=0;q2=0;q3=0;

end 
 
quz=sqrt(q1^2+q2^2+q3^2);
n1=A/sqrt(A^2+B^2+C^2);
n2=B/sqrt(A^2+B^2+C^2);
n3=C/sqrt(A^2+B^2+C^2);
kap=(q1*A+B*q2+C*q3)/sqrt(A^2+B^2+C^2);
d= kap^2*(A^2+B^2+C^2)/(a^2*A^2+b^2*B^2+c^2*C^2);

ak=1;
bk=-(n1^2*(1/b^2+1/c^2)+n2^2*(1/c^2+1/a^2)+n3^2*(1/b^2+1/a^2));
ck=(n1/b/c)^2+(n2/a/c)^2+(n3/b/a)^2;
p=[ak bk ck];
kok=roots(p);
Bye=sqrt((1-d)/kok(1));
Aye=sqrt((1-d)/kok(2));


%%fprintf('Ellipse of intersection half length axis\n a=%12.4f  b=%12.4f  \n center coordinates of ellips q1=%12.4f q2=%12.4f q3=%12.4f\n',Aye,Bye,q1,q2,q3)
Ayy=Aye;Byy=Bye;
if length(elp)==3
Tmat=eye(4);
end
ym=Tmat*[q1 q2 q3 1]';qx=ym(1);qy=ym(2);qz=ym(3);
end
function [param,artik] = duzdenkbelirleme( XYZ )
%
format long
[m,n]=size(XYZ);
x=-inv(XYZ)*ones(m,1);
param=x;
artik=XYZ*x+ones(m,1);
end

function [center,radii,quat,R] = ellipsoid_im2ex(v)
% Cast ellipsoid defined with implicit parameter vector to explicit form.
% The implicit equation of a general ellipse is
% F(x,y,z) = Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz - 1 = 0
%
% Input arguments:
% v:
%    the 10 parameters describing the ellipsoid algebraically
% Output arguments:
% center:
%    ellispoid center coordinates [cx; cy; cz]
% ax:
%    ellipsoid semi-axes (radii) [a; b; c]
% quat:
%    ellipsoid rotation in quaternion representation
% R:
%    ellipsoid rotation (radii directions as rows of the 3x3 matrix)
%
% See also: ellipse_im2ex

% Copyright 2011 Levente Hunyadi

validateattributes(v, {'numeric'}, {'nonempty','real','vector'});

% eliminate times two from rotation and translation terms
v = v(:);
validateattributes(v, {'numeric'}, {'size',[10,1]});
v(4:9) = 0.5*v(4:9);

% find the algebraic form of the ellipsoid (quadratic form matrix)
Q = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];

% find the center of the ellipsoid
center = Q(1:3,1:3) \ -v(7:9);

if nargout > 1
    % form the corresponding translation matrix
    T = eye(4,4);
    T(4, 1:3) = center;

    % translate to the center
    S = T * Q * T';

    % check for positive definiteness
    [~,indef] = chol( -S(4,4)*S(1:3,1:3) );
    if indef > 0  % matrix is not positive definite
        error('ellipsoid_im2ex:InvalidArgumentValue', ...
            'Parameters do not define a real ellipse.');
    end
    
    % solve the eigenproblem
    [evecs, evals] = eig( S(1:3,1:3) );
    radii = realsqrt( -S(4,4) ./ diag(evals) );
    %radii = sqrt( -S(4,4) ./ diag(evals) ); %ben deðiþtirdim

    % convert rotation matrix to quaternion
    if nargout > 2
        quat = rot2quat(evecs);
        R = evecs';
    end
end
end