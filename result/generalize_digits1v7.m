digit_bestV = load('digits1v7_bestV.mat');
digit_bestV = digit_bestV.v_best;
digit_bestB = load('digits1v7_bestB.mat');
digit_bestB = digit_bestB.b_best;

sample = load('optdigits.tes');
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



C = 0.5;
d = 1;
sigma = 5;

sample = load('optdigits.tra');
train = zeros(0,n);
for i=1:m
    if sample(i,n+1)==1 || sample(i,n+1)==7
        train(end+1,:) = sample(i,1:n);
    end
end
cluster = zeros(m,1);
for i=1:m
    train(end,:) = data(i,:);
    linear_kernel = train * train';  % 线性核函数
    polynomial_kernel = (train * train'+1).^d;  % 多项式核函数
    guass_kernel = zeros(m,m);
    for i = 1:m1
        for j = 1:m1
            data_temp = train(i,:)-train(j,:);
            guass_kernel(i,j) = exp(-data_temp*data_temp'/(2*sigma));  % 高斯核函数
        end
    end
    X = {};
    X{1} = getFeature(linear_kernel);
    X{2} = getFeature(polynomial_kernel);
    X{3} = getFeature(guass_kernel);

    temp = zeros(1,m);
    for i=1:3
        X_temp = X{i};
        temp = temp + digit_bestV(:,i)' * X_temp(end,:)';
    end
    temp = temp + digit_bestB;
    cluster(end+1) = sign(temp);
end
