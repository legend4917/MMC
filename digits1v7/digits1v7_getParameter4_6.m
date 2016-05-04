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

C_all = [];    % �ɳڱ������Ʋ���
for i=-3:-1
    C_all(end+1) = 2^i;
end
sigma_all = [0.5 1 2 5 7 10 12 15 17 20];   % ��˹�˲���
d_all = [];     % ����ʽ�˲���
for i=0:4
    d_all(end+1) = 2^i;
end
l = m/3;  % ��ƽ��Լ������
M = 3;  % ʹ�õĺ˺�������
alpha = 0.01;   % SOCP������ֹ����
epsilon = 0.05;
acc = zeros(200,5);


num = 1;
linear_kernel = data * data';  % ���Ժ˺���
for C_num=1:length(C_all)   % �������ŵ��ɳڱ������Ʋ���
    C = C_all(C_num);
    for d_num=1:length(d_all)   % �������ŵĶ���ʽ�˲��� 
        d = d_all(d_num);
        polynomial_kernel = (data * data'+1).^d;  % ����ʽ�˺���
        for sigma_num=1:length(sigma_all)   % �������ŵĸ�˹�˲���
            sigma = sigma_all(sigma_num);
            guass_kernel = zeros(m,m);
            for i = 1:m
                for j = 1:m
                    data_temp = data(i,:)-data(j,:);
                    guass_kernel(i,j) = exp(-data_temp*data_temp'/(2*sigma));  % ��˹�˺���
                end
            end
            X = {};
            X{1} = getFeature(linear_kernel);
            X{2} = getFeature(polynomial_kernel);
            X{3} = getFeature(guass_kernel);
            c_arr = {};         % ����������Υ����Լ��
            acc_arr = [];
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
                            -l <= sum(l_temp+b) <= l;    % ����Լ��
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
                fprintf('��ǰ�����(%d)����Υ����Լ��\n',numel(c_arr));
                temp = zeros(1,m);
                for i=1:M
                    temp = temp + v(:,i)' * X{i}';
                end
                temp = (temp + b)';
                c = zeros(m,1);     % ��ʼ����Υ��Լ����Ӧ��c����
                xi_m = 0;
                for i = 1:m     % Ѱ����Υ����Լ��
                    if abs(temp(i)) < 1
                        c(i) = 1;
                    else
                        c(i) = 0;
                    end
                    xi_m = xi_m + c(i) * (1 - abs(temp(i)));
                end
                xi_m = xi_m / m;
                fprintf('��ƽ��(%f)\n',xi_m-xi);
                fprintf('����(%f, %f, %f)\n', C_num, d_num, sigma_num);
                if xi_m-xi <= epsilon
                    cluster = sign(temp);
                    break;
                end
                c_arr{numel(c_arr)+1} = c;
                cluster = sign(temp);
                acc_arr(end+1) = sum(label==cluster)/m;
                if acc_arr(end) < 0.5
                    acc_arr(end) = 1 - acc_arr(end);
                end
                disp(acc_arr);
                if length(acc_arr) >= 150
                    break;
                end
            end
            acc(num,1) = C;
            acc(num,2) = d;
            acc(num,3) = sigma;
            acc_arr(end+1) = sum(label==cluster)/m;
            if acc_arr(end) < 0.5
                acc_arr(end) = 1 - acc_arr(end);
            end
            acc(num,4) = max(acc_arr);
            disp(['��ƽ����������������Ϊ��',num2str(num)]);
            num = num + 1;
        end
    end
end