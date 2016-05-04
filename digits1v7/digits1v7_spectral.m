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


similarityFunction = inline('exp(-sum((vecA-vecB).^2)/(2*sigma^2))','vecA','vecB','sigma');

similarityMat_origin = zeros(m,m);
for i=1:m
    for j=i:m
        similarityMat_origin(i,j) = similarityFunction(data(i,:),data(j,:),0.5);
        similarityMat_origin(j,i) = similarityMat_origin(i,j);
    end
end

acc = 0;
for iter = 1:50
    acc_spectral = [];
    epsilon = 10^(-100);
    for ep=1:100
        epsilon = epsilon * 10;
        similarityMat = zeros(m,m);
        for i=1:m
            for j=i:m
                if similarityMat_origin(i,j) < epsilon
                    similarityMat(i,j) = 0;
                else
                    similarityMat(i,j) = similarityMat_origin(i,j);
                end
                similarityMat(j,i) = similarityMat(i,j);
            end
        end

        D = zeros(m,m);
        for i=1:m
            D(i,i) = sum(similarityMat(i,:));
        end

        L = D - similarityMat;
        [P,lambda] = eig(L);
        data = P(:,1:2); 
        [idx,C] = kmeans(data,2);
        idx = (idx - 1) * 2 - 1;

        acc_spectral_temp = sum(label==idx(:,1))/m;
        if acc_spectral_temp < 0.5
            acc_spectral_temp = 1 - acc_spectral_temp;
        end
        acc_spectral(end+1) = acc_spectral_temp;
    end
    acc = acc + max(acc_spectral);
    disp(max(acc_spectral))
end

acc = acc / 50;