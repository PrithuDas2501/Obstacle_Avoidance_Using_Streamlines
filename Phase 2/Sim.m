clear all; clc;

WP_ObstacleData;
%%
% model = multirotor;
% s = state(model);
% s(end) = 0;
Time_Of_Sim = 0:0.05:2;
[t,X] = ode45(@(t,X) dState(t, X, pp, A, c), [0 tim], [WP(:,1);0;0;0]);
coord = ppval(pp,t);
%%
K = X;
% %%
% figure(1); clf;
% plot(t,K(:,3))
% figure(2); clf;
% plot(t,K(:,2))
% figure(3); clf;
% plot(t,K(:,1))
% %% Pitch
% figure(4); clf;
% plot(t,K(:,8))
% %% Roll
% figure(5); clf;
% plot(t,K(:,9))
%%
figure(6);
plot3(K(:,1),K(:,2),K(:,3),'g', LineWidth=2)
hold on
xlabel('X');
ylabel('Y');
zlabel('Z');
r = svd(A);
[X,Y,Z] = ellipsoid(0,0,0,r(1)^-0.5,r(2)^-0.5,r(3)^-0.5);
s = surf(c(1)+X,c(2)+Y,c(3)+Z,FaceColor='b',FaceAlpha=0.1);
scatter3(WP(1,:), WP(2,:), WP(3,:),100, 'k','filled')
plot3(coord(1,:),coord(2,:),coord(3,:),'r',LineWidth=2)
grid on

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
%%
gammap1 = zeros(length(K),1);
for i = 1:length(K)
    gammap1(i) = findgamma(A,c,[K(i,1);K(i,2);K(i,3)]);
end
figure(7);
hold on
grid on
xlabel('Time')
ylabel('\Gamma')
plot(t,gammap1,'g',LineWidth=1);
%%
%scatter3(0,1,0,60,'b')
% %%
% figure(7); clf;
% hold on
% grid on
% legend()
% plot(t,180*K(:,10)/pi,LineWidth=1.5,Color = 'g', DisplayName='Roll Rate (p)')
% plot(t,180*K(:,11)/pi,LineWidth=1.5,Color = 'k', DisplayName='Pitch Rate (q)')
% %plot(t,K(:,12),LineWidth=1.5,DisplayName='Yaw Rate (r)')
% xlabel('Time (s)')
% ylabel('Body Axis Turn Rates (degree/s)')
% title('Turn Rates')
% %%
% figure(8); clf;
% hold on
% grid on
% legend()
% plot(t,180*K(:,8)/pi,LineWidth=1.5,Color='k', DisplayName='Pitch Angle')
% plot(t,180*K(:,9)/pi,LineWidth=1.5,Color='g', DisplayName='Roll Angle')
% %plot(t,K(:,10),LineWidth=1.5,DisplayName='Yaw Angle')
% xlabel('Time (s)')
% ylabel('Rotation Angles (degree)')
% title('Turn Angles')
