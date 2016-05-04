clear;
clc;
sample = load('optdigits.tra');
[m,n] = size(sample);
data = zeros(0,n-1);
label = [];
for i=1:m
    if sample(i,n)==1
        data(end+1,:) = sample(i,1:n-1);
        label(end+1) = 1; 
    end
    if sample(i,n)==7
        data(end+1,:) = sample(i,1:n-1);
        label(end+1) = -1; 
    end
end
label = label';
[m,n] = size(data);


acc_kmean = 0;
for i=1:50
    [idx,C] = kmeans(data,2);

    idx = (idx - 1) * 2 - 1;
    temp = sum(idx==label)/m;
    if temp < 0.5
        temp = 1 - temp;
    end
    acc_kmean = acc_kmean + temp;
end
acc_kmean = acc_kmean / 50;

