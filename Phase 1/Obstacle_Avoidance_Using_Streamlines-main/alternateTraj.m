function y = alternateTraj(Traj,streamval,alpha, rx, ry, Zdash, Xdashdash, Ydashdash,Pos)
%% Find Projection
proj = Traj - (dot(Zdash,Traj) - dot(Zdash,Pos))*Zdash;
options = optimset('Display','off');
Coord = fsolve(@(X) Inplane(X,Xdashdash,Ydashdash,proj), [0;0],options);
ycoord = fsolve(@(Y) gety(Y,Coord(1),streamval,alpha,rx,ry),0,options);
xcoord = Coord(1);

y = xcoord*Xdashdash + ycoord*Ydashdash;
end