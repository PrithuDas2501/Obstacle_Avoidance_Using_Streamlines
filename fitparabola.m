% Loosely based on  matlabcentral answers 80541
%https://www.mathworks.com/matlabcentral/answers/80541
% 
%  This algorithm does a RANSAC-type search over possible rotation angles
%  to find the angle at which the sum-square residual fit is minimized.
%  Because noisy data often ends up finding a minimum which is essentially
%  a straight line along the parabola's axis of symmetry, a second test is
%  run at 90 degrees to the first fit.  Whichever of these two yields the
%  larger coefficient of x^2 is deemed the correct fit. A final
%  minimization is done in a local neighborhood to optimize the fit.
%
% This can be sensitive to large amounts of noise in the input, and often
% fails for noisy data when only one "leg" of the parabola is involved. 
% One thing that can help is to observe the the source data. If it
% looks like a  horizontal parabola (x = y^2) rotate90 deg BEFORE
%  running this fitting function.
%
%   INPUTS
%       x: vector of x-coords, OR a 2xN or Nx2 matrix.  If x is a matrix,
%       then the y-input is ignored.
%       y:  y-coordinates of the dataset
%       thetaguess: (not required.) angle in radians. Provides a
%        starting point for the    optimization routine. 
%
%   OUTPUT
%       pfit: a structure containing
%           abc -- 3 values as ax^2 + bx + c = y for the parabola when
%           rotated to vertical
%           vertex -- vertex of the parabola when rotated to vertical
%           theta -- the angle, radians the source is rotated from vertical
%           resid -- the norm of the residuals at best fit.
%       
%   SEE ALSO: 
%       fitellipse
%       fithyperbola
%
function pfit  = fitparabola(x,y,thetaguess)
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

xy=[x(:),y(:)].';
if ~exist('thetaguess','var'),
    thetaguess = 0 ; 
end
% input validation
thetaguess = thetaguess(1); 
coefftmp = zeros(3,1); % 'cause ONLY doing 2nd order polynom
        thetatmp = fminsearch(@(theta) pcost(theta,xy), thetaguess );    
        [costtmp,coefftmp(:,1)]=pcost(thetatmp,xy);

     mincost = min(costtmp);
     maxcost = max(costtmp);
%      thewinner = find(costtmp == mincost);
%      thewinner = thewinner(1) ; % just in case multiples
     coeffs = coefftmp;%   
     theta = thetatmp;  %  
     resid = [mincost,maxcost]; % to get an idea of quality 
% after finding cost min,  Just try theta+90
theta90 = theta - pi/2;
[cost90,coeffs90] = pcost(theta90,xy);
   [a90,b90,c90]=deal(coeffs90(1),coeffs90(2), coeffs90(3));
    abc90 = [a90,b90,c90]; % I want to know everything! 
        xv90 = -b90/2/a90; % this is h in vertex form, k = f(h) .
     vertex90=R(-theta90)*[xv90;polyval(abc90,xv90)];

 % a,b,c are the ax^2 +bx+c coeffs unrotated
% deal() just distributes the coeffs.. same as [a,b,c] = coeffs{:} 
    [a,b,c]=deal(coeffs(1),coeffs(2), coeffs(3));
    abc = [a,b,c]; % I want to know everything! 
        xv = -b/2/a; % this is h in vertex form, k = f(h) .
  % 'vertex form' is y = a(x-h)^2 + k and [h,k] is vertex coords
 %polyval  evaluates 2nd arg at polynom coeffs set  
     vertex=R(-theta)*[xv;polyval(coeffs,xv)];
% various diagnostics currently commented out
 pfit.vertex = vertex;
 pfit.theta = theta;
 pfit.abc  = abc;
 pfit.resid = resid;
%  pfit.vertex90 = vertex90;
%  pfit.theta90 = theta90;
%  pfit.abc90  = abc90;
%  pfit.resid90 = cost90;
% 
%  secondary optimization
% first, compare the "a" coeffs at theta and theta90. the larger magnitude
% is the one we want.  If it's at theta, we're done. If not, optimize near
% theta90 .  Don't think it's possible to have a third local minimum.
if abs(a90) > abs(a),
    theta2 = fminbnd(@(theta) pcost(theta,xy),theta90-pi/4,...
        theta90+pi/4 );    
    [pfit.resid, pfit.abc] = pcost(theta2, xy);
    pfit.theta = theta2;
    xv2 = -pfit.abc(2)/2/pfit.abc(1); 
    pfit.vertex = R(-theta2)*[xv2;polyval(pfit.abc,xv2)];
end
    
end

% Enable these functions when posting to FEX
 function [Cost,coeffs,Rxy] = pcost(theta,xy)
   Rxy=R(theta) * xy;
 % Rxy = rotcoord(theta,xy);
   [coeffs,S] = polyfit(Rxy(1,:),Rxy(2,:),2);
    Cost=S.normr; %norm of the residuals. R calls this lm$residuals
 end
 
 function Rmat=R(theta)
     Rmat=[cos(theta), -sin(theta); sin(theta), cos(theta)];
 end
	 
     