function y = alternateTraj(Interpx,Interpy,d,xbarbar,ybarbar,t,T)
if t>=T
    point2d = [Interpx(t);Interpy(t)];
    y = d + [xbarbar, ybarbar]*point2d;
else
    point2d = [Interpx(T);Interpy(T)];
    y = d + [xbarbar, ybarbar]*point2d;
end
end