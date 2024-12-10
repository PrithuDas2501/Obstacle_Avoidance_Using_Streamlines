StartPoint = [0; 1; 0];

WP = [0  10 50 60;
      0  0 0 5;
      0  5 10 10];

% WP = [0 0 10 10;
%       0 0 0 10;
%       0 10 10 10];
%WP = [0 25;
%      1 25;
%      0 25];
WP_Tolerance = 1;


%% Rectangle Obstacle
h = 10;
l = 6;
Obs_Loc = [30; 0; 0];
Obs_P_Matrix = [Obs_Loc(1)-l/2 Obs_Loc(1)-l/2 Obs_Loc(1)+l/2 Obs_Loc(1)+l/2 Obs_Loc(1)-l/2 Obs_Loc(1)-l/2 Obs_Loc(1)+l/2 Obs_Loc(1)+l/2;
                Obs_Loc(2)-l/2 Obs_Loc(2)+l/2 Obs_Loc(2)+l/2 Obs_Loc(2)-l/2 Obs_Loc(2)-l/2 Obs_Loc(2)+l/2 Obs_Loc(2)+l/2 Obs_Loc(2)-l/2;
                0              0              0              0              h              h              h              h             ];
%%
P = Obs_P_Matrix;
Tol = 0.1;
[A,c] = MinVolEllipse(P,Tol);
% A = [1/4 0 0; 0 1/9 0; 0 0 1/25]/100;
%%
tim = 200;
[q,qd,qdd,qddd,qdddd,pp,tPoints,tSamples] = minsnappolytraj(WP,0:tim/(size(WP,2)-1):tim,10000);
dpp = ppDer(pp);       % Velocity
ddpp = ppDer(dpp);     % Acceleration
dddpp = ppDer(ddpp);   % Jerk
ddddpp = ppDer(dddpp); % Snap