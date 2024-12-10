% Adapted from a paper on fitting ellipse, hyperbola
% Journal of Electronic Imaging / July 2004 / Vol. 13(3) / 495
%  Paul O'Leary, Paul Zsombor-Murray

%
%   INPUTS
%       x: vector of x-coords, OR a 2xN or Nx2 matrix.  If x is a matrix,
%       then the y-input is ignored.
%       y: vector of y-coords
%
%   OUTPUTS
%       efit contains the following:
%       -- a vector of six coefficients representing, in order,
%       Ax^2 + Bxy + Cy^2 + Dx + Ey + F.   Note that these can be converted
%       into [center, axes, rotation angle] via the  AtoG.m function.
%       At this time, fit quality parameters such as rms or rss error are
%       not provided; these are easily calculated from the coefficient set
%
%   SEE ALSO: 
%       fitparabola
%       fithyperbola
%

function efit  =  fitellipse(x, y)

% some notes from original author:
%  "need to make data zero-mean and scale before fitting "
% Prior to starting the fit procedure, a translation is applied to
% place the centroid of the data at the origin, then isotropic
% scaling is made to ensure that the root mean square distance
% of the points to the origin is &. The matrix corresponding
% to this similarity transform is:
% S = [s,0, -s<x>; 0,s, -s<y>; 0,0,1] 
% where s = n*sqrt(2)/ (sum(i=1:n) sqrt(x(i)-<x>)^2+(y(i)-<y>)^2)
% so feed (input*S) to func, and apply S^-1 *answer to descale
%
if size(x,2) == 2, % vals in columns
    y = x(:,2);
    x = x(:,1);
elseif size(x,1) == 2, % vals in rows
    y = x(2,:)';
    x = x(1,:)';
else
    x = x(:);
    y = y(:); % guarantee columns
end
if length(x) ~= length(y)
    error('x and y length mismatch')
end
x2  =  x.^2;
y2  =  y.^2;
xy  =  x.*y;
%%Set up the design and scatter matrices
%
D1 = [x2,xy,y2];
D2 = [x,y,ones(size(x))];
%
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
%%test the rank of S3
%
[Us3, Ss3, Vs3] = svd( S3 ); % 
condNrs = diag( Ss3 )/Ss3(1,1); % note this is scalar division, not matrix
%
%epsilon = 1e–10;
epsilon = eps;
if condNrs(3)<epsilon
warning('S3 is degenerate');
return;
end;
%define the constraint matrix and its inverse
C = [0, 0, -2;
0, 1, 0;
-2, 0, 0];
% Ci = inv(C);
%Setup and solve the generalized eigenvector problem
%
T = -inv( S3 )*S2';
% mathworks says left-div is better than inv * mat
S = C\(S1 + S2 * T);
% S = Ci * (S1 + S2 * T);
%
[evec,eval] = eig( S );
%%evaluate and sort resulting constraint values
%
cond = evec(2,:).^2 - 4*(evec(1,:).*evec(3,:));
[condVals index] = sort( cond );
%
% NB  alpha1,alpha2 = a,b,c,d,e,f for standard conic equation. (x^2,xy,y^2,x,y,1)
% Pick up the elliptical solution
%
eValE = condVals( 1 );
alpha1 = evec( :,index(1) );
alpha2 = T*alpha1;
efit = [alpha1;alpha2];
%%Pick up the hyperbolic solution
%
% possibleHs = condVals(2:3)+condVals(1);
% [minDiffVal, minDiffAt] = min( abs( possibleHs ) );
% eValH = possibleHs(minDiffAt);
% alpha1 = evec( :,index(minDiffAt+1) );
% alpha2 = T*alpha1;
% hyperbola = [alpha1; alpha2];
end
