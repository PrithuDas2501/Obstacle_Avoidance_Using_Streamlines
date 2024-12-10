function y = TrajfromM(M,Axe,Bye)
y = [];
[a,l] = size(M);
i = 1;
rem = l;
while rem>0
    partition = M(2,i);
    check = checkinsidellipse(Axe,Bye,M(1,i+1),M(2,i+1));
    if check>0
    y = [y M(:,i+1:i + partition)];
    end
    rem = rem - partition - 1;
    i = i+1+partition;
end
end