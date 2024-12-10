clear; clc;
%%
model = multirotor;
s0 = state(model);
% s0(6) = -1;
e = environment(model);
lim = 20*pi/180;
Wt = model.Configuration.Mass*e.Gravity;
model.Configuration.PDRoll(2) = 0;
model.Configuration.PDPitch(2) = 0;
WP_ObstacleData;
%%
tic;
[t,S] = ode45(@(t,s) DSTATE(t, s, model, e, Wt, pp, dpp, ddpp, dddpp, ddddpp, A, c, lim), [0 tim], zeros(12,1));
simtime = toc
timescale = simtime/tim
%%
coord = ppval(pp,t);
TrajVel = ppval(dpp,t);
TrajAccel = ppval(ddpp,t);
TrajJerk = ppval(dddpp,t);
% TrajSnap = ppval(ddddpp,t);
%%
figure(1); clf
hold on
grid on
xlabel('t')
ylabel('X Error')
plot(t,coord(1,:)'-S(:,1))
%%
figure(2); clf
hold on
grid on
xlabel('t')
ylabel('Y Error')
plot(t,coord(2,:)'-S(:,2))
%%
figure(3); clf
hold on
grid on
xlabel('t')
ylabel('Z Error')
plot(t,coord(3,:)'-S(:,3))
%%
figure(4); clf
hold on
grid on
legend('on')
xlabel('t')
ylabel('Thrust')
K_dz = 3.5;
K_pz = 3.8;
Thrust = Wt + 0.47*(K_dz*(TrajVel(3,:)'-S(:,6)) + K_pz*(coord(3,:)'-S(:,3)));
for i = 1:length(Thrust)
    if Thrust(i)<0
        Thrust(i)=0;
    end
end
plot(t,Thrust,DisplayName='Thrust Commanded');
plot([0 tim], [Wt Wt], 'r:',DisplayName='Weight');
%%
figure(5); clf;
hold on
grid on
legend('on')
xlabel('t')
ylabel('Euler Angles')
plot(t, S(:,7),DisplayName='\psi')
plot(t, S(:,8),DisplayName='\theta')
plot(t, S(:,9),DisplayName='\phi')
plot([0 tim], [-lim -lim], 'r:',DisplayName='Limit');
plot([0 tim], [lim lim], 'r:',DisplayName='Limit');
%%
figure(6); clf;
hold on
grid on
legend('on')
xlabel('t')
ylabel('Rotation Rates')
plot(t, S(:,10),DisplayName='p')
plot(t, S(:,11),DisplayName='q')
plot(t, S(:,12),DisplayName='r')
plot([0 tim], [-lim -lim], 'r:',DisplayName='Limit');
plot([0 tim], [lim lim], 'r:',DisplayName='Limit');
%%%%%%%%%%%%%%%%%% 
ControlEff = trapz(t,S(:,10).^2+S(:,11).^2+S(:,12).^2)
%%
figure(7); clf
hold on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3(S(:,1), S(:,2), S(:,3),'g',LineWidth=2)
r = svd(A);
[X,Y,Z] = ellipsoid(0,0,0,r(1)^-0.5,r(2)^-0.5,r(3)^-0.5);
s = surf(c(1)+X,c(2)+Y,c(3)+Z,FaceColor='b',FaceAlpha=0.1);
scatter3(WP(1,:), WP(2,:), WP(3,:),100, 'k','filled')
plot3(coord(1,:),coord(2,:),coord(3,:),'r',LineWidth=2)

%Bottom Surface
[X,Y] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),Obs_P_Matrix(2,1):Obs_P_Matrix(2,2));
Z = 0*X;
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
%Top Surface
[X,Y] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),Obs_P_Matrix(2,1):Obs_P_Matrix(2,2));
Z = 0*X + h;
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
%Side Surfaces
[X,Z] = meshgrid(Obs_P_Matrix(1,1):0.5:Obs_P_Matrix(1,3),0:h);
Y = 0*Z + Obs_P_Matrix(2,1);
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
Y = 0*Z + Obs_P_Matrix(2,2);
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
[Y,Z] = meshgrid(Obs_P_Matrix(2,1):0.5:Obs_P_Matrix(2,2),0:h);
X = 0*Z + Obs_P_Matrix(1,1);
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
X = 0*Z + Obs_P_Matrix(1,3);
surface(X,Y,Z,FaceColor='k',FaceAlpha=0.7)
% axis equal
%%
gammap1 = zeros(length(S),1);
for i = 1:length(S)
    gammap1(i) = findgamma(A,c,[S(i,1);S(i,2);S(i,3)]);
end
figure(8); clf
hold on
grid on
xlabel('Time')
ylabel('\Gamma')
plot(t,gammap1,'g',LineWidth=1); 
%%
function y = DSTATE(t, s, model, e, Wt, pp, dpp, ddpp, dddpp, ddddpp, A, c, lim)
% s(3) = -s(3);
% s(6) = -s(6);
persistent CollisionFlag Axe Bye d xbarbar ybarbar stream2d Interpx Interpy Interpvx Interpvy T GAMMAprev accelvec k0_Interp Tol;
if isempty(CollisionFlag)
    CollisionFlag = 0;
    T = 0;
    Tol = 7;
    k0_Interp = EllipticalGainSchedule(Tol,0.8*Tol);
    GAMMAprev = findgamma(A,c,[s(1);s(2);s(3)]);
    accelvec=[0;0;-9.81];
end

if CollisionFlag==0
    T=t;
end

GAMMA = findgamma(A,c,[s(1);s(2);s(3)]); 

Distance = distancePointLine3D([s(1);s(2);s(3)],c,[0;0;1]);
if GAMMA <Tol
    options = optimoptions('fmincon', 'Display', 'off');
    [lambda,fval] = fmincon(@(lambda) CheckCollosion(lambda, A, c, [s(1);s(2);s(3)], [s(4);s(5);s(6)]),0, [], [],[],[],[],[],[],options);

    if lambda >0 && fval<0 && CollisionFlag == 0
        CollisionFlag = 1;
        DesiredSlope1 = atand(-s(4)/(sign(s(5))*sqrt(max(0.01^2,s(5)^2))));
        VOL = [cosd(DesiredSlope1); sind(DesiredSlope1); 0];
        VOR = -VOL;
        VelVec = [s(4);s(5);s(6)]/max(0.001,norm([s(4);s(5);s(6)]));
        VOV = cross(VOR,VelVec);
        VerticalCoeff = 0;
        normal = VerticalCoeff*VOL + (1-VerticalCoeff)*VOV;
        plane = [normal' -normal'*[s(1);s(2);s(3)]];
        % normal = fmincon(@(normal) OptimPlane(normal,x,c,a),normal);
        % normal = normal/norm(normal);
        % plane = [normal' -normal'*x(1:3)'];
        options = optimoptions('fmincon','display','none');
        optimset('display','off');
        d = fmincon(@(d) find_d(d,c,A),c,[],[],normal',-plane(4));

        qmin = fmincon(@(q) findmin(q,d,c,A),[s(1);s(2);s(3)],[],[],normal',-plane(4));
        kmin = fsolve(@(k) ((d + k*(qmin-d) - c)'*A*(d + k*(qmin-d) - c)-1)^2,1,options);
        qmin = d + kmin*(qmin-d);
        Bye = norm(qmin-d);

        qmaj = d + cross(qmin-d,normal);
        kmaj = fsolve(@(k) ((d + k*(qmaj-d) - c)'*A*(d + k*(qmaj-d) - c)-1)^2,1,options);
        qmaj = d + kmaj*(qmaj-d);
        Axe = norm(qmaj-d);

        xbarbar = (qmaj-d)/norm(qmaj-d);
        ybarbar = (qmin-d)/norm(qmin-d);

        xbases = [xbarbar, ybarbar]\([s(1);s(2);s(3)]-d);
        Velbases = [xbarbar, ybarbar]\[s(4);s(5);s(6)];

        Nontransformedtheta = -atan((Velbases(2)/Velbases(1))*Axe/Bye);
        Nontransformedxbases = [Bye/Axe 0;0 1]*xbases;
        NontransformedVel = [Bye/Axe 0;0 1]*Velbases;

        rotmat = [cos(Nontransformedtheta) -sin(Nontransformedtheta);
                  sin(Nontransformedtheta)  cos(Nontransformedtheta)];

        NontransformedNonrotatedxbases = rotmat*Nontransformedxbases;
        NontransformedNonrotatedVel = rotmat*NontransformedVel;
        z1 = (NontransformedNonrotatedxbases(1,1)/Bye + 1i*NontransformedNonrotatedxbases(2,1)/Bye);
        psi = imag(z1+1./z1);
        xlim = NontransformedNonrotatedxbases(1,1);
        ylim = sign(NontransformedNonrotatedxbases(2,1))*max(2*Bye,abs(NontransformedNonrotatedxbases(2,1)));
        Altx = linspace(xlim,-xlim,500);
        Alty = linspace(ylim,-ylim,500);
        [xgrid, ygrid] = meshgrid(Altx,Alty);
        z = (xgrid/Bye + 1i*ygrid/Bye);
        psifunc = imag(z+1./z);
        M=contourc(Altx,Alty,psifunc,[psi,psi]);
        %M0=contourc(Altx,Alty,psifunc,[sign(psi)*0.01,sign(psi)*0.01]);
        stream2d = TrajfromM_(M,Axe,Bye,NontransformedNonrotatedxbases,NontransformedNonrotatedVel,Nontransformedtheta); %[xs;ys;px;py;vx;vy;t]
        % Interpx = griddedInterpolant(stream2d(7,:)+t, stream2d(1,:));
        % Interpy = griddedInterpolant(stream2d(7,:)+t, stream2d(2,:));
        % Interpvx = griddedInterpolant(stream2d(7,:)+t, stream2d(3,:));
        % Interpvy = griddedInterpolant(stream2d(7,:)+t, stream2d(4,:));
        Interpx = griddedInterpolant(stream2d(7,:)+t, stream2d(8,:));
        Interpy = griddedInterpolant(stream2d(7,:)+t, stream2d(9,:));
        Interpvx = griddedInterpolant(stream2d(7,:)+t, stream2d(10,:));
        Interpvy = griddedInterpolant(stream2d(7,:)+t, stream2d(11,:));
        disp(t);
        %stream2d_0 = TrajfromM_(M0,Axe,Bye,NontransformedNonrotatedxbases,NontransformedNonrotatedVel,Nontransformedtheta); %[xs;ys;px;py;vx;vy;t]
    end
% else
%     if CollisionFlag==1
%         if Distance == 5
%             CollisionFlag=0;
%         end
%     end
end

%% This this is giving issues, why is t reducing with iterations??
if CollisionFlag==1
    if t-T>stream2d(7,end)
        CollisionFlag=0;
    end
end

if CollisionFlag==1
    %% Insert alternate Traj
    P1 = alternateTraj(Interpx,Interpy,d,xbarbar,ybarbar,t,T);
    P0 = ppval(pp,t);
    P1dot = alternateTraj(Interpvx,Interpvy,[0;0;0],xbarbar,ybarbar,t,T);
    P0dot = ppval(dpp,t);

    % if GAMMA>GAMMAprev
    %     k0 = (2*sigmoid(GAMMA^1)-1);
    % else
    %     k0 = 0;
    % end

    k0 = (2*sigmoid(abs(GAMMA^0.8))-1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % k0 = cos_based_k0(GAMMA,Tol);
    % k0 = GAMMA/Tol;
    % if k0>1
    %     k0=1;
    % end

    %k0 = k0_Interp(GAMMA);

    Traj = k0*P0 + (1-k0)*P1;
    TrajVel = k0*P0dot + (1-k0)*P1dot;
    TrajAccel = ppval(ddpp,t);
    TrajJerk = ppval(dddpp,t);
    % TrajAccel = [0;0;0];
    % TrajJerk = [0;0;0];
    % 
    % Traj = alternateTraj(Interpx,Interpy,d,xbarbar,ybarbar,t,T);
    % TrajVel = alternateTraj(Interpvx,Interpvy,d,xbarbar,ybarbar,t,T);
    % TrajAccel = ppval(ddpp,t);
    % TrajJerk = ppval(dddpp,t);
else
    Traj = ppval(pp,t);
    TrajVel = ppval(dpp,t);
    TrajAccel = ppval(ddpp,t);
    TrajJerk = ppval(dddpp,t);
end

GAMMAprev = GAMMA;
% t
% Traj = ppval(pp,t);
% TrajVel = ppval(dpp,t);
% TrajAccel = ppval(ddpp,t);
% TrajJerk = ppval(dddpp,t);
% % TrajSnap = ppval(ddddpp,t);

%% Quadcopter Model
% Constants
K_dz = 3.5;
K_pz = 3.8;
K_dphi = 0.2;% 0.5;
K_pphi = 2;
K_dtheta = 0.2;% 0.2;
K_ptheta = 2;% 1.8;
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

% Desired Values
xdotdotc = TrajAccel(1) + K_px*(Traj(1)-s(1)) + K_dx*(TrajVel(1)-s(4));
ydotdotc = TrajAccel(2) + K_py*(Traj(2)-s(2)) + K_dy*(TrajVel(2)-s(5));

xdotdotdotc = TrajJerk(1) + K_pjx*(TrajVel(1)-s(4)) + K_djx*(TrajAccel(1)-accelvec(1)); %%
ydotdotdotc = TrajJerk(2) + K_pjy*(TrajVel(2)-s(5)) + K_djy*(TrajAccel(2)-accelvec(2)); %%

Yawc = 0;
Rollc = (xdotdotc*sin(Yawc)-ydotdotc*cos(Yawc))/9.81;
Pitchc = (xdotdotc*cos(Yawc)+ydotdotc*sin(Yawc))/9.81;

if Pitchc > lim
    Pitchc = lim;
elseif Pitchc < -lim
    Pitchc = -lim;
end
if Rollc > lim
   Rollc = lim;
elseif Rollc < -lim
    Rollc = -lim;
end

Yawdotc = 0;
Rolldotc = (xdotdotdotc*sin(Yawc)+xdotdotc*cos(Yawc)*Yawdotc-ydotdotdotc*cos(Yawc)+ydotdotc*sin(Yawc)*Yawdotc)/9.81;
Pitchdotc= (xdotdotdotc*cos(Yawc)-xdotdotc*sin(Yawc)*Yawdotc+ydotdotdotc*sin(Yawc)+ydotdotc*cos(Yawc)*Yawdotc)/9.81;

bodyratesc = [1  0         -sin(s(8));
              0  cos(s(9))  sin(s(9))*cos(s(8));
              0 -sin(s(9))  cos(s(9))*cos(s(8))]*[Rolldotc;Pitchdotc;Yawdotc]; % [pc;qc;rc]
% 
% bodyratesc = eye(3)*[Rolldotc;Pitchdotc;Yawdotc];

% Control Input
u1 = Wt + 0.47*(K_dz*(TrajVel(3)-s(6)) + K_pz*(Traj(3)-s(3)));
if u1<0
    u1=0;
end
u2 = K_dphi*(bodyratesc(1)-s(10)) + K_pphi*(Rollc-s(9));
u3 = K_dtheta*(bodyratesc(2)-s(11)) + K_ptheta*(Pitchc-s(8));
u4 = K_dpsi*(bodyratesc(3)-s(12)) + K_ppsi*(Yawc-s(7));
U = [u1;u2;u3;u4];

% State Variables
y = zeros(12,1);
y(1) = s(4);
y(2) = s(5);
y(3) = s(6);

RotMat = [cos(s(8))*cos(s(7))   cos(s(7))*sin(s(9))*sin(s(8))-cos(s(9))*sin(s(7))    cos(s(9))*cos(s(7))*sin(s(8))+sin(s(9))*sin(s(7));
          cos(s(8))*sin(s(7))   cos(s(9))*cos(s(7))+sin(s(7))*sin(s(8))*sin(s(9))   -cos(s(7))*sin(s(9))+cos(s(9))*sin(s(8))*sin(s(7));
         -sin(s(8))             cos(s(8))*sin(s(9))                                  cos(s(9))*cos(s(8))];

accelvec = ([0;0;-Wt] + RotMat*[0;0;U(1)])/0.47;
y(4) = accelvec(1);
y(5) = accelvec(2);
y(6) = accelvec(3);

y(7) = s(12);
y(8) = s(11);
y(9) = s(10);

J_p = 0.0086;
J_q = 0.0086;
J_r = 0.0176;

y(10) = U(2)/J_p;
y(11) = U(3)/J_q;
y(12) = U(4)/J_r;
end