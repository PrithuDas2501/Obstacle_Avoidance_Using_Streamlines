function y = QuadcopterModel(t, s, Wt, pp, dpp, ddpp, dddpp, ddddpp)
%% Trajectory Elements
Traj = ppval(pp,t);
TrajVel = ppval(dpp,t);
TrajAccel = ppval(ddpp,t);
TrajJerk = ppval(dddpp,t);
% TrajSnap = ppval(ddddpp,t);

%% Constants
K_dz = 3.5;
K_pz = 3.8;
K_dphi = 0.5;
K_pphi = 12.8;
K_dtheta = 0.2;
K_ptheta = 1.8;
K_dpsi = 0.5;
K_ppsi = 2;
K_djx = 3.5;
K_pjx = 6;
K_djy = 4.2;
K_pjy = 12.7;
K_dx = 3.5;
K_px = 6;
K_dy = 4.2;
K_py = 12.7;

%% Desired Values
xdotdotc = TrajAccel(1) + K_px*(Traj(1)-s(1)) + K_dx*(TrajVel(1)-s(4));
ydotdotc = TrajAccel(2) + K_py*(Traj(2)-s(2)) + K_dy*(TrajVel(2)-s(5));

xdotdotdotc = TrajJerk(1) + K_pjx*(TrajVel(1)-s(4)) + K_djx*(TrajAccel(1)-accelvec(1)); %%
ydotdotdotc = TrajJerk(2) + K_pjy*(TrajVel(2)-s(5)) + K_djy*(TrajAccel(2)-accelvec(2)); %%

Yawc = 0;
Rollc = (xdotdotc*sin(Yawc)-ydotdotc*cos(Yawc))/9.81;
Pitchc = (xdotdotc*cos(Yawc)+ydotdotc*sin(Yawc))/9.81;

Yawdotc = 0;
Rolldotc = (xdotdotdotc*sin(Yawc)+xdotdotc*cos(Yawc)*Yawdotc-ydotdotdotc*cos(Yawc)+ydotdotc*sin(Yawc)*Yawdotc)/9.81;
Pitchdotc= (xdotdotdotc*cos(Yawc)-xdotdotc*sin(Yawc)*Yawdotc+ydotdotdotc*sin(Yawc)+ydotdotc*cos(Yawc)*Yawdotc)/9.81;

bodyratesc = [1  0         -sin(s(8));
              0  cos(s(9))  sin(s(9))*cos(s(8));
              0 -sin(s(9))  cos(s(9))*cos(s(8))]*[Rolldotc;Pitchdotc;Yawdotc]; % [pc;qc;rc]

%% Control Input
u1 = Wt + 0.47*(K_dz*(TrajVel(3)-s(6)) + K_pz*(Traj(3)-s(3)));
u2 = K_dphi*(bodyratesc(1)-s(10)) + K_pphi*(Rollc-s(9));
u3 = K_dtheta*(bodyratesc(2)-s(11)) + K_ptheta*(Pitchc-s(8));
u4 = K_dpsi*(bodyratesc(3)-s(12)) + K_ppsi*(Yawc-s(7));
U = [u1;u2;u3;u4];

%% State Variables
% y = zeros(13,1);
y = zeros(12,1);
y(1) = s(4);
y(2) = s(5);
y(3) = s(6);

RotMat = [cos(s(8))*cos(s(7))   cos(s(7))*sin(s(9))*sin(s(8))-cos(s(9))*sin(s(7))    cos(s(9))*cos(s(7))*sin(s(8))+sin(s(9))*sin(s(7));
          cos(s(8))*sin(s(7))   cos(s(9))*cos(s(7))+sin(s(7))*sin(s(8))*sin(s(9))   -cos(s(7))*sin(s(9))+cos(s(9))*sin(s(8))*sin(s(7));
         -sin(s(8))             cos(s(8))*sin(s(9))                                  cos(s(9))*cos(s(8))];

accelvec = ([0;0;Wt] + RotMat*[0;0;U(1)])/0.47;
y(4) = accelvec(1);
y(5) = accelvec(2);
y(6) = accelvec(3);

% J = [1 sin(s(9))*tan(s(8)) cos(s(9))*tan(s(8));
%      0 cos(s(9))          -sin(s(9));
%      0 sin(s(9))/cos(s(8)) cos(s(9))/cos(s(8))];
% eulrate = J*[s(10);s(11);s(12)];
% 
% y(7) = eulrate(1);
% y(8) = eulrate(2);
% y(9) = eulrate(3);

y(7) = s(10);
y(8) = s(11);
y(9) = s(12);

J_p = 0.047;
J_q = 0.047;
J_r = 0.08;

y(10) = U(2)/J_p;
y(11) = U(3)/J_q;
y(12) = U(4)/J_r;
end