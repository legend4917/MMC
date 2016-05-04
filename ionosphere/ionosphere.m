clear;
clc;
sample = load('ionosphere.data.txt');
[m,n] = size(sample);   % 样本个数与维度
data = sample(:,1:n-1);
label = sample(:,n);
[m,n] = size(data);

C = 0.0625;
d = 4;
sigma = 5;
epsilon = 0.01;



l = 105;  % 类平衡约束条件
M = 3;  % 使用的核函数个数
alpha = 0.01;   % SOCP问题终止条件
error = 0.1;

linear_kernel = data * data';  % 线性核函数
polynomial_kernel = (data * data'+1).^d;  % 多项式核函数
guass_kernel = zeros(m,m);
for i = 1:m
    for j = 1:m
        data_temp = data(i,:)-data(j,:);
        guass_kernel(i,j) = exp(-data_temp*data_temp'/(2*sigma));  % 高斯核函数
    end
end
X = {};
X{1} = getFeature(linear_kernel);
X{2} = getFeature(polynomial_kernel);
X{3} = getFeature(guass_kernel);
c_arr = {};         % 保存每轮CCCP收敛后的最违背约束
acc_arr = [];       % 保存每轮CCCP收敛后的精度
xi_arr = [];        % 保存每轮CCCP收敛后的xi值
time_arr = [];
t1 = clock;
while 1
    v_k = ones(m,M);
    b_k = 1;
    res_temp = inf;
    while 1
        cvx_begin sdp quiet
            cvx_solver mosek
            cvx_save_prefs
            variables betta(M) T(M)
            variables xi(1) b(1)
            variable v(m,M)
            minimize(0.5*sum(T)+C*xi)
            subject to
                for i=1:M
                    norm([2*v(:,i);T(i)-betta(i)]) <= T(i)+betta(i);
                    betta(i) >= 0;
                end
                sum(betta.^2) <= 1;
                xi >= 0;
                l_temp = zeros(1,m);
                for i=1:M
                    l_temp = l_temp + v(:,i)'*X{i}';
                end
                -l <= sum(l_temp+b) <= l;    % 公共约束
                disp(numel(c_arr));
                if numel(c_arr) > 0
                    z_temp = zeros(1,m);
                    for i=1:M
                        z_temp = z_temp + v_k(:,i)'*X{i}';
                    end
                    z = sign(z_temp+b_k)';
                    for i=1:numel(c_arr)
                        1/m*sum(c_arr{i}) <= xi+1/m*((c_arr{i}.*z)'*(l_temp+b)');
                    end
                end
        cvx_end
        res = cvx_optval;
        if abs((res-res_temp)/res) <= alpha
            break;
        end
        v_k = v;
        b_k = b;
        res_temp = res;
    end
    temp = zeros(1,m);
    for i=1:M
        temp = temp + v(:,i)' * X{i}';
    end
    temp = (temp + b)';
    c = zeros(m,1);     % 初始化最违背约束对应的c向量
    xi_m = 0;
    for i = 1:m     % 寻找最违背的约束
        if abs(temp(i)) < 1
            c(i) = 1;
        else
            c(i) = 0;
        end
        xi_m = xi_m + c(i) * (1 - abs(temp(i)));
    end
    xi_m = xi_m / m;
    xi_arr(end+1) = xi_m - xi;
    while xi_m - xi <= error
        t2 = clock;
        time_arr(end+1) = etime(t2,t1);
        error = error - 0.01;
    end
    fprintf('割平面(%f)\n',xi_m-xi);
    cluster = sign(temp);
    acc_arr(end+1) = sum(label==cluster)/m;
    if acc_arr(end) < 0.5
        acc_arr(end) = 1 - acc_arr(end);
    end
    disp(max(acc_arr));
    if xi_m-xi <= epsilon
        break;
    end
    c_arr{numel(c_arr)+1} = c;
end
acc = max(acc_arr);