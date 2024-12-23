clear%, close
%% Initializing Figure Window
figure(1); 
hold on
grid on 
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal
options = optimset('Display','off');
%% Ellipsoid
rx = 3;
ry = 2;
rz = 5;
c = [0;0;0];
a = diag([1/rx^2,1/ry^2,1/rz^2]);
eul = [0 0 0];
rotm = eul2rotm(eul);
eulZYX = rotm2eul(rotm);
quat = eul2quat(eulZYX);
[X,Y,Z] = ellipsoid(0,0,0,rx,ry,rz);
s = surf(c(1)+X,c(2)+Y,c(3)+Z,FaceColor='r',FaceAlpha=0.1);
if acosd(quat(1)) ==0
else
    rotate(s,quat(2:4),2*acosd(quat(1)))
end
%% Cutting Planes
%x = [5,4,7,-5,-5,-5];
x = [6,6,0,-5,-4,1];
VertCoeff = 0;i=1;
%VertCoeff = linspace(0,1,4);
%for i = 1:4
    
    DesiredSlope1 = atand(-x(4)/(sign(x(5))*sqrt(max(0.01^2,x(5)^2))));
    VOL = [cosd(DesiredSlope1); sind(DesiredSlope1); 0];
    VOR = -VOL;
    VelVec = [x(4);x(5);x(6)]/max(0.001,norm([x(4);x(5);x(6)]));
    VOV = cross(VOR,VelVec);
    VerticalCoeff = VertCoeff(i);
    normal = VerticalCoeff*VOL + (1-VerticalCoeff)*VOV;
    plane = [normal' -normal'*x(1:3)'];
    % normal = fmincon(@(normal) OptimPlane(normal,x,c,a),normal);
    % normal = normal/norm(normal);
    % plane = [normal' -normal'*x(1:3)'];
    options = optimoptions('fmincon','display','none');
    optimset('display','off');
    d = fmincon(@(d) find_d(d,c,a),c,[],[],normal',-plane(4));
    
    qmin = fmincon(@(q) findmin(q,d,c,a),x(1:3)',[],[],normal',-plane(4));
    kmin = fsolve(@(k) ((d + k*(qmin-d) - c)'*a*(d + k*(qmin-d) - c)-1)^2,1,options);
    qmin = d + kmin*(qmin-d);
    Bye = norm(qmin-d);
    
    qmaj = d + cross(qmin-d,normal);
    kmaj = fsolve(@(k) ((d + k*(qmaj-d) - c)'*a*(d + k*(qmaj-d) - c)-1)^2,1,options);
    qmaj = d + kmaj*(qmaj-d);
    Axe = norm(qmaj-d);
    
    xbarbar = (qmaj-d)/norm(qmaj-d);
    ybarbar = (qmin-d)/norm(qmin-d);
    
    scatter3(d(1),d(2),d(3),50,"yellow","filled");
    scatter3(x(1),x(2),x(3),50,"black","filled");
    scatter3(qmaj(1),qmaj(2),qmaj(3),50,"Green","filled");
    scatter3(qmin(1),qmin(2),qmin(3),50,"Blue","filled");
    
    xlim = linspace(c(1) - 10*VelVec(1),c(1) + 10*VelVec(1),2);
    ylim = linspace(c(2) - 10*VelVec(2),c(2) + 10*VelVec(2),2);
    p1 = d+xlim(1)*xbarbar + ylim(1)*ybarbar;
    p2 = d+xlim(2)*xbarbar + ylim(1)*ybarbar;
    p3 = d+xlim(2)*xbarbar + ylim(2)*ybarbar;
    p4 = d+xlim(1)*xbarbar + ylim(2)*ybarbar;
    P = [p1,p2,p3,p4];
    fill3(P(1,:),P(2,:),P(3,:),'b', 'FaceAlpha', 0.1)
%%
xbases = [xbarbar, ybarbar]\(x(1:3)'-d);
Velbases = [xbarbar, ybarbar]\VelVec(1:3);
% Nontransformedtheta = -atan((Velbases(2)/Velbases(1))*Axe/Bye);
% e = cos(Nontransformedtheta) + 1i*sin(Nontransformedtheta);
% z1 = (xbases(1,1)/Axe + 1i*xbases(2,1)/Bye)*e;
% psi = imag(z1+1./z1)
% % Altx = linspace(xbases(1,1),-xbases(1,1),1000);
% % Alty = linspace(xbases(2,1),-xbases(2,1),1000);
% xlim = sign(xbases(1,1))*max(2*Axe,abs(xbases(1,1)));
% ylim = sign(xbases(2,1))*max(2*Bye,abs(xbases(2,1)));
% Altx = linspace(xlim,-xlim,1000);
% Alty = linspace(ylim,-ylim,1000);
% [xgrid, ygrid] = meshgrid(Altx,Alty);
% z = (xgrid/Axe + 1i*ygrid/Bye)*e;
% psifunc = imag(z+1./z);
% M=contourc(Altx,Alty,psifunc,[psi,psi]);
% stream2d = TrajfromM(M,Axe,Bye);
% figure(2)
% hold on
% grid on 
% axis equal
% l = linspace(0,2*pi,500);
% plot(Axe*cos(l),Bye*sin(l))
% %plot(M(1,2:1+M(2,1)), M(2,2:1+M(2,1)),'k',LineWidth=1);
% plot(stream2d(1,:),stream2d(2,:),'k',LineWidth=1);
%%
Nontransformedtheta = -atan((Velbases(2)/Velbases(1))*Axe/Bye);
Nontransformedxbases = [Bye/Axe 0;0 1]*xbases;
NontransformedVel = [Bye/Axe 0;0 1]*Velbases;
e = cos(Nontransformedtheta) + 1i*sin(Nontransformedtheta);
% z1 = (Nontransformedxbases(1,1) + 1i*Nontransformedxbases(2,1))*e;
% psi = imag(z1+1./z1);

rotmat = [cos(Nontransformedtheta) -sin(Nontransformedtheta);
          sin(Nontransformedtheta)  cos(Nontransformedtheta)];

NontransformedNonrotatedxbases = rotmat*Nontransformedxbases;
NontransformedNonrotatedVel = rotmat*NontransformedVel;
z1 = (NontransformedNonrotatedxbases(1,1)/Bye + 1i*NontransformedNonrotatedxbases(2,1)/Bye);
psi = imag(z1+1./z1);
% Altx = linspace(xbases(1,1),-xbases(1,1),1000);
% Alty = linspace(xbases(2,1),-xbases(2,1),1000);
% xlim = sign(xbases(1,1))*max(2*Axe,abs(xbases(1,1)));
% ylim = sign(xbases(2,1))*max(2*Bye,abs(xbases(2,1)));
% xlim = xbases(1,1);
% ylim = xbases(2,1);
xlim = NontransformedNonrotatedxbases(1,1);
ylim = sign(NontransformedNonrotatedxbases(2,1))*max(2*Bye,abs(NontransformedNonrotatedxbases(2,1)));
Altx = linspace(xlim,-xlim,500);
Alty = linspace(ylim,-ylim,500);
[xgrid, ygrid] = meshgrid(Altx,Alty);
z = (xgrid/Bye + 1i*ygrid/Bye);
psifunc = imag(z+1./z);
M=contourc(Altx,Alty,psifunc,[psi,psi]);
M0=contourc(Altx,Alty,psifunc,[sign(psi)*0.01,sign(psi)*0.01]);
stream2d = TrajfromM_(M,Axe,Bye,NontransformedNonrotatedxbases,NontransformedNonrotatedVel,Nontransformedtheta); %[xs;ys;px;py;vx;vy;t]
%stream2d_0 = TrajfromM_(M0,Axe,Bye,NontransformedNonrotatedxbases,NontransformedNonrotatedVel,Nontransformedtheta); %[xs;ys;px;py;vx;vy;t]
%%
figure(2); clf
hold on
grid on 
axis equal
l = linspace(0,2*pi,500);
plot(Axe*cos(l),Bye*sin(l))
%plot(M(1,2:1+M(2,1)), M(2,2:1+M(2,1)),'k',LineWidth=1);
plot(stream2d(1,:),stream2d(2,:),'r',LineWidth=1);
plot(stream2d(8,:),stream2d(9,:),'g',LineWidth=1);
%%
figure(1); hold on
%stream2d = M(:,2:1+M(2,1));
%stream2d = M;
stream3d = d + [xbarbar, ybarbar]*stream2d(1:2,:);
stream3d_better = d + [xbarbar, ybarbar]*stream2d(8:9,:);
ellipse = d + [xbarbar, ybarbar]*[Axe*cos(l); Bye*sin(l)];
plot3(stream3d(1,:),stream3d(2,:),stream3d(3,:),'r',LineWidth=2);
plot3(stream3d_better(1,:),stream3d_better(2,:),stream3d_better(3,:),'g',LineWidth=2);
plot3(ellipse(1,:),ellipse(2,:),ellipse(3,:),'k',LineWidth=2);
%end