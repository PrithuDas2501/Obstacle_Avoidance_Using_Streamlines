function out = TrajfromM_(M,Axe,Bye,xbases,Vel,Nontransformedtheta)
y = [];
[a,l] = size(M);
i = 1;
rem = l;
while rem>0
    partition = M(2,i);
    check = checkinsidellipse(Bye,Bye,M(1,i+1),M(2,i+1));
    if check>0
    y = [y M(:,i+1:i + partition)];
    end
    rem = rem - partition - 1;
    i = i+1+partition;
end

% if norm(y(:,1)-xbases)>0.5
%     y = fliplr(y);
% end
% len = length(y);
% minnorm = 100;
% for j = 1:len
%     if j == 1
%         minnorm = (xbases-y(:,j))'*(xbases-y(:,j));
%         continue
%     else
%         n = (xbases-y(:,j))'*(xbases-y(:,j));
%         if n<=minnorm
%             minnorm = n;
%         else
%             break
%         end
%     end
% end
% 
% y = y(:,j:end);

% Proj = [Vel(2,1)/Vel(1,1) -1;Vel(1,1)/Vel(2,1)
% 1]\[(xbases(1,1)*Vel(2,1)/Vel(1,1) - xbases(2,1))*ones(1,length(y));
% y(1,:)*Vel(1,1)/Vel(2,1) + y(2,:)]; Singularity issues with this one
Proj = [y(1,:);y(2,1)*ones(1,length(y))];
num = length(Proj);
Tseq = (Proj - xbases*ones(1,num))./Vel; Tseq = Tseq(1,:);

len = length(y);
minnorm = 100;
for j = 1:len
    if j == 1
        minnorm = (xbases-Proj(:,j))'*(xbases-Proj(:,j));
        continue
    else
        n = (xbases-Proj(:,j))'*(xbases-Proj(:,j));
        if n<=minnorm
            minnorm = n;
        else
            break
        end
    end
end

y = y(:,j:end);
Tseq = Tseq(:,j:end);
Proj = Proj(:,j:end);
num = length(Tseq);
V = Vel*ones(1,num);
V_ = Vel*ones(1,num);
V_2= Vel*ones(1,num);
d1 = abs(-Vel(2)*xbases(1)/Vel(1)+xbases(2))/sqrt(1+(Vel(2)/Vel(1))^2);
d2=sqrt(y(1,ceil(num/2))^2+y(2,ceil(num/2))^2);
scale = (Bye-d1)/(d2-d1);
scaled_y = zeros(2,num);
scaled_y2= zeros(2,num);
for i = 1:num
    theta = atand(-Vel(1)/Vel(2));
    scaled_y(:,i) = Proj(:,i) - scale*(Vel(2)*(y(1,i)-xbases(1))/Vel(1)-y(2,i) + xbases(2))/sqrt(1+(Vel(2)/Vel(1))^2)*[cosd(theta);sind(theta)];
    scaled_y2(:,i) = Proj(:,i) + scale*(Vel(2)*(y(1,i)-xbases(1))/Vel(1)-y(2,i) + xbases(2))/sqrt(1+(Vel(2)/Vel(1))^2)*[cosd(theta);sind(theta)];
    if i>1
        V(:,i) = (y(:,i)-y(:,i-1))./(Tseq(i)-Tseq(i-1));
        V_(:,i) = (scaled_y(:,i)-scaled_y(:,i-1))./(Tseq(i)-Tseq(i-1));
        V_2(:,i) = (scaled_y2(:,i)-scaled_y2(:,i-1))./(Tseq(i)-Tseq(i-1));
    end
end

rotmat = [Axe/Bye 0;0 1]*[cos(Nontransformedtheta) sin(Nontransformedtheta);-sin(Nontransformedtheta)  cos(Nontransformedtheta)];
%rotmat = [Axe/Bye 0;0 1];
y = rotmat*y;
Proj = rotmat*Proj;
V = rotmat*V;
scaled_y = rotmat*scaled_y;
V_ = rotmat*V_;
scaled_y2 = rotmat*scaled_y2;
V_2 = rotmat*V_2;

if norm(scaled_y-y)<norm(scaled_y2-y)
    out = [y;Proj;V;Tseq;scaled_y;V_];
else
    out = [y;Proj;V;Tseq;scaled_y2;V_2];
end

end