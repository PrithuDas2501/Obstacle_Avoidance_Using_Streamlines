function y = OptimPlane(theta,VOL,VOV,x,c,a)
options = optimoptions('fmincon','display','none');
normal = cos(theta)*VOL + sin(theta)*VOV;
plane = [normal' -normal'*x(1:3)'];
d = fmincon(@(d) find_d(d,c,a),c,[],[],normal',-plane(4));

qmin = fmincon(@(q) findmin(q,d,c,a),x(1:3)',[],[],normal',-plane(4));
kmin = fsolve(@(k) ((d + k*(qmin-d) - c)'*a*(d + k*(qmin-d) - c)-1)^2,1,options);
qmin = d + kmin*(qmin-d);
Bye = norm(qmin-d);

qmaj = d + cross(qmin-d,normal);
kmaj = fsolve(@(k) ((d + k*(qmaj-d) - c)'*a*(d + k*(qmaj-d) - c)-1)^2,1,options);
qmaj = d + kmaj*(qmaj-d);
Axe = norm(qmaj-d);

y = Axe+Bye;
end