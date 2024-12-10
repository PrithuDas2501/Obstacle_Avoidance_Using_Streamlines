x = linspace(0,10,1000);
y = zeros(1,1000);
shift = 4;
for i = 1:1000
    % y(i) = sigmoid(x(i)-shift);
    % y(i) = 0.5*(1-cos(x(i)));
    y(i) = (2*sigmoid(abs(x(i)^1))-1);
end
figure(10)
hold on
grid on
xlabel('\Gamma')
ylabel('K_0')
title('K_0 Graph')
plot(x,y,'g')
% plot([shift shift], [-1 1], 'r:');
% scatter(shift,0,40,'red')