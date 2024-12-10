function y = dState(t, X, pp, A, c)


persistent CollisionFlag Axe Bye d xbarbar ybarbar stream2d Interpx Interpy T GAMMAprev;
if isempty(CollisionFlag)
    CollisionFlag = 0;
    T = 0;
    GAMMAprev = findgamma(A,c,[X(1);X(2);X(3)]);
    %figure;
    %grid on;
    %hold on;
end

if CollisionFlag==0
    T=t;
end

GAMMA = findgamma(A,c,[X(1);X(2);X(3)]); 

Distance = distancePointLine3D(X(1:3),c,[0;0;1]);
if Distance <12
    options = optimoptions('fmincon', 'Display', 'off');
    [lambda,fval] = fmincon(@(lambda) CheckCollosion(lambda, A, c, [X(1);X(2);X(3)], [X(4);X(5);X(6)]),0, [], [],[],[],[],[],[],options);

    if lambda >0 && fval<0 && CollisionFlag == 0
        CollisionFlag = 1;
        DesiredSlope1 = atand(-X(4)/(sign(X(5))*sqrt(max(0.01^2,X(5)^2))));
        VOL = [cosd(DesiredSlope1); sind(DesiredSlope1); 0];
        VOR = -VOL;
        VelVec = [X(4);X(5);X(6)]/max(0.001,norm([X(4);X(5);X(6)]));
        VOV = cross(VOR,VelVec);
        VerticalCoeff = 1;
        normal = VerticalCoeff*VOL + (1-VerticalCoeff)*VOV;
        plane = [normal' -normal'*X(1:3)];
        % normal = fmincon(@(normal) OptimPlane(normal,x,c,a),normal);
        % normal = normal/norm(normal);
        % plane = [normal' -normal'*x(1:3)'];
        options = optimoptions('fmincon','display','none');
        optimset('display','off');
        d = fmincon(@(d) find_d(d,c,A),c,[],[],normal',-plane(4));

        qmin = fmincon(@(q) findmin(q,d,c,A),X(1:3),[],[],normal',-plane(4));
        kmin = fsolve(@(k) ((d + k*(qmin-d) - c)'*A*(d + k*(qmin-d) - c)-1)^2,1,options);
        qmin = d + kmin*(qmin-d);
        Bye = norm(qmin-d);

        qmaj = d + cross(qmin-d,normal);
        kmaj = fsolve(@(k) ((d + k*(qmaj-d) - c)'*A*(d + k*(qmaj-d) - c)-1)^2,1,options);
        qmaj = d + kmaj*(qmaj-d);
        Axe = norm(qmaj-d);

        xbarbar = (qmaj-d)/norm(qmaj-d);
        ybarbar = (qmin-d)/norm(qmin-d);

        xbases = [xbarbar, ybarbar]\(X(1:3)-d);
        Velbases = [xbarbar, ybarbar]\[X(4);X(5);X(6)];

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
        Interpx = griddedInterpolant(stream2d(7,:)+t, stream2d(8,:));
        Interpy = griddedInterpolant(stream2d(7,:)+t, stream2d(9,:));
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
    if GAMMA>GAMMAprev
        k0 = (2*sigmoid(GAMMA^1)-1);
    else
        k0 = 0;
    end
    Traj = k0*P0 + (1-k0)*P1;
    % Traj = alternateTraj(Interpx,Interpy,d,xbarbar,ybarbar,t,T);
else
    Traj = ppval(pp,t);
end

GAMMAprev = GAMMA;
%plot3(Traj(1),Traj(2),Traj(3));
% u = model.control;
% en = environment(model);
% 
% s = zeros(13,1);
% s(1:12) = X(1:end-1);
% 
% u.Thrust = model.Configuration.Mass*en.Gravity/(cos(s(8))*(cos(s(9))));
% u.Thrust = u.Thrust + (1*(Traj(3)+X(3)) + 0.1*(0+X(6)))/(cos(s(8))*(cos(s(9))));
% % desYaw = atan((Traj(2)-X(2))/(Traj(1)-X(1)));
% % if anynan(desYaw)==1
% %     desYaw=pi/4;
% % end
% % u.YawRate = -1*(X(7)-desYaw);
% % relx = [cos(desYaw) sin(desYaw);-sin(desYaw) cos(desYaw)]*[Traj(1)-X(1);Traj(2)-X(2)];
% % relv = [cos(desYaw) sin(desYaw);-sin(desYaw) cos(desYaw)]*[-X(4);-X(5)];
% % u.Pitch = u.Pitch - 1*relx(1) - 1*relv(1);
% u.Roll = u.Roll + 1*(Traj(2)-X(2)) + 0.2*(0-X(5));
% u.Pitch = u.Pitch + 1*(Traj(1)-X(1)) + 0.2*(0-X(4));
% 
% if u.Thrust>1.2*model.Configuration.Mass*en.Gravity
%     u.Thrust = 1.2*model.Configuration.Mass*en.Gravity;
% end
% lim = pi/18;
% if u.Pitch > lim
%     u.Pitch = lim;
% elseif u.Pitch < -lim
%     u.Pitch = -lim;
% end
% if u.Roll > lim
%     u.Roll = lim;
% elseif u.Roll < -lim
%     u.Roll = -lim;
% end
% 
% s(end) = u.Thrust;
% y = derivative(model,s,u,en);

y(1,1) = X(4);
y(2,1) = X(5);
y(3,1) = X(6);

y(4,1) = (0.5*(Traj(1)-X(1)) + 1*(0-X(4)));
y(5,1) = (0.5*(Traj(2)-X(2)) + 1*(0-X(5)));
y(6,1) = (0.5*(Traj(3)-X(3)) + 1*(0-X(6)));
t
end