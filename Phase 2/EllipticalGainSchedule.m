function y = EllipticalGainSchedule(Gamma0,Gamma_Half)

    n = 1000;
    
    x1 = linspace(Gamma_Half,Gamma0,n);
    y1 = zeros(1,n);
    x2 = linspace(0,Gamma_Half,n);
    y2 = zeros(1,n);
    
    for i = 1:n
        y1(i) = elli(x1(i),Gamma0,0.5,Gamma0-Gamma_Half,0.5);
        y2(i) = elli(x2(i),0,0.5,Gamma_Half,0.5);
    end
    y = griddedInterpolant([x2(1:end-1),x1], [1-y2(1:end-1),y1]);

end

% xplot = linspace(0,Gamma0,100);
% yplot = zeros(1,100);
% for i = 1:100
%     yplot(i) = Out(xplot(i));
% end
% %%
% figure
% hold on
% plot(xplot, yplot,'b',LineWidth=2)
% % plot(x1,y1)
% % plot(x2,1-y2)