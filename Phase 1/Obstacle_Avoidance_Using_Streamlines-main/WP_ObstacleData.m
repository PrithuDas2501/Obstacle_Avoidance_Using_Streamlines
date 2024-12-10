StartPoint = [0; 1; 0];

WP = [0  0   6   18  0;
      1  10  20  0   0;
      0 -20 -20 -20 -20];
WP_Tolerance = 1;

%% Rectangle Obstacle
h = 40;
l = 3;
Obs_Loc = [14; 14; 0];
Obs_P_Matrix = [Obs_Loc(1)-l/2 Obs_Loc(1)-l/2 Obs_Loc(1)+l/2 Obs_Loc(1)+l/2 Obs_Loc(1)-l/2 Obs_Loc(1)-l/2 Obs_Loc(1)+l/2 Obs_Loc(1)+l/2;
                Obs_Loc(2)-l/2 Obs_Loc(2)+l/2 Obs_Loc(2)+l/2 Obs_Loc(2)-l/2 Obs_Loc(2)-l/2 Obs_Loc(2)+l/2 Obs_Loc(2)+l/2 Obs_Loc(2)-l/2;
                0              0              0              0              h              h              h              h             ];
%%
P = Obs_P_Matrix;
Tol = 0.1;
[A,c] = MinVolEllipse(P,Tol);
%A = [1/4 0 0; 0 1/9 0; 0 0 1/25]/100;
%%
[q,qd,qdd,qddd,qdddd,pp,tPoints,tSamples] = minsnappolytraj(WP,0:20:80,1000);
