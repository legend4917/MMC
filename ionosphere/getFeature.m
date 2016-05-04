function [ X ] = getFeature( K )
n = length(K);
[p,x] = eig(K);
for i=1:n
    if x(i,i) < 0
        x(i,i) = 0;
    end
    x(i,i) = x(i,i)^0.5;
end
X = p*x;
end