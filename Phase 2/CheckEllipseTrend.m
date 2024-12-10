n = 50;
theta = linspace(0,pi,n);
Array = zeros(1,n);
for i = 1:n
    Array(1,i) = OptimPlane(theta(i),VOL,VOV,x,c,a);
end
%%
figure
hold on
grid on
xlabel('Normal Direction')
ylabel('A + B')
title('Sum of axes length across different plane intersections')
plot(theta*180/pi,Array(1,:))
