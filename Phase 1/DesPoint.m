function vec = DesPoint(t)
T = ones(10,1);
T = t*T;
for i = 1:10
    T(i) = t(i)^(10-i);
end

DesPoint = 
end