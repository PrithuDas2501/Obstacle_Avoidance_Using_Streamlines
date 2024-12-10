function y = FindAOA(theta, Rx, C, Xdash, Zdash, A)
%% Needs to account for non centre aligned velocity!!!!!
%% Actual Angle of attack will be this angle found minus angle between
Rotation_Matrix = RotMatrix(theta,Zdash);
y = ((C+Rx*Rotation_Matrix*Xdash)'*A*(C+Rx*Rotation_Matrix*Xdash)-1)^2;
end