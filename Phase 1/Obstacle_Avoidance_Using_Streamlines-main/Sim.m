clear ; clc;
WP_ObstacleData;
%%
model = multirotor;
s = state(model);
s(end) = 0;
Time_Of_Sim = 0:0.05:2;
[t,X] = ode23(@(t,X) dState(t, X, model, pp, A, c), [0 60], s);
%%
K = X;
%%
figure
plot(t,K(:,3))
figure
plot(t,K(:,2))
figure
plot(t,K(:,1))
%% Pitch
figure
plot(t,K(:,8))
%% Roll
figure
plot(t,K(:,9))
%%
figure
plot3(K(:,1),K(:,2),K(:,3))
hold on
xlabel('X');
ylabel('Y');
zlabel('Z');
r = svd(A);
%ellipsoid(c(1),c(2),c(3),r(1),r(2),r(3))
grid on
%%
%Bottom Surface
[X,Y] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),Obs_P_Matrix(2,1):Obs_P_Matrix(2,2));
Z = 0*X;
surface(X,Y,Z)
%Top Surface
[X,Y] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),Obs_P_Matrix(2,1):Obs_P_Matrix(2,2));
Z = 0*X + h;
surface(X,Y,Z)
%Side Surfaces
[X,Z] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),0:h);
Y = 0*Z + Obs_P_Matrix(2,1);
surface(X,Y,Z)
Y = 0*Z + Obs_P_Matrix(2,2);
surface(X,Y,Z)
[Y,Z] = meshgrid(Obs_P_Matrix(2,1):0.5:Obs_P_Matrix(2,2),0:h);
X = 0*Z + Obs_P_Matrix(1,1);
surface(X,Y,Z)
X = 0*Z + Obs_P_Matrix(1,3);
surface(X,Y,Z)
