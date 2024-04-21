% 要求日均总能量最大，只需以年为周期，求年最大

% 如果使用fmincon求最大值
f = @(inclines) -calculate_year_energy(inclines(1), inclines(2)); % 注意这里取负值(因为要求最大值，后续会取相反数)
x0 = [20, 30]; % 初始猜测值
A = []; % 线性不等式约束矩阵
b = []; % 线性不等式约束向量
Aeq = []; % 线性等式约束矩阵
beq = []; % 线性等式约束向量
lb = [0, -90]; % 变量下界
ub = [90, 90]; % 变量上界
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_max1,f_valmax1] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, [], options);

% 如果结合遗传算法求最大值
fitnessFunction = @(inclines) -calculate_year_energy(inclines(1), inclines(2));
numberOfVariables = 2; % 变量的数量
% 遗传算法选项设置
options = optimoptions('ga', 'MutationFcn', {@mutationuniform, 0.05}); % 设置突变函数和突变率，这里的0.05是突变率
options = optimoptions(options, 'Display', 'iter', 'PopulationSize', 100, 'MaxGenerations', 100, 'CrossoverFraction', 0.8); % 其他参数设置
% 执行遗传算法
[x_max2, f_valmax2] = ga(fitnessFunction, numberOfVariables, A, b, Aeq, beq, lb, ub, [], options);

% 计算平均能量值，由于fitnessFunction返回负值，所以这里取反
s = 1; % 光伏板面积
f_val_average1 = -f_valmax1 / 365 * s;
f_val_average2 = -f_valmax2 / 365 * s;

% 输出结果
disp('fmincon求得最大值点：');
disp(x_max1);
disp('对应的平均日能量值：');
disp(f_val_average1);
disp('结合遗传算法求得最大值点：');
disp(x_max2);
disp('对应的平均日能量值：');
disp(f_val_average2);

% 绘制二元函数图像
% 定义x和y的范围
x = linspace(23, 26, 15); % 生成一个线性间隔的向量
y = linspace(-1.5, 1.5, 15);
% 创建网格坐标矩阵
[X, Y] = meshgrid(x, y);
% 计算每个点上的函数值
Z = zeros(15,15);
for i=1:15
    for j=1:15
        Z(i,j) = calculate_year_energy(x(i), y(j));
    end
end
surf(X, Y, Z) % 绘制三维曲面图
title('每平方米太阳辐射量W/m^2与水平仰角和方位角的关系') % 添加标题和轴标签
xlabel('水平仰角')
ylabel('方位角')
zlabel('每平方米太阳辐射量W/m^2')
grid on % 开启网格
colorbar % 添加颜色条

% 计算一年单位面积光伏板吸收的辐射值
function energy_of_year = calculate_year_energy(incline_1, incline_2)
    monthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
    energy_1 = [1405 1394 1378 1353 1334 1316 1308 1315 1330 1350 1372 1392];
    energy_of_year = 0;
    k = 0.000401466772612056;
    for i=1:12
        for j=1:monthdays(i)
            f1 = @(t1)calculate_height(2025, i, j, t1, 30.58333333333);
            % 计算日出日落时间
            t_min = fzero(f1, [3, 8.5]);
            t_max = fzero(f1, [16, 21]);
            f2 = @(t2)calculate_energy_3(2025, i, j, t2, incline_1, incline_2, k, energy_1(i));
            energy_of_year = energy_of_year + 3600 * integral(f2, t_min, t_max);
        end
    end
end

% 计算太阳高度角
function hangle = calculate_height(year,month,day,time,latitude)
monthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
latitude = deg2rad(latitude); % 转化纬度为弧度制
day_of_year = day;
for i = 1:month-1
    day_of_year = day_of_year+monthdays(i); % 计算年积日
end
N0 = 79.6764+0.2422*(year-1985)-floor((year-1985)/4); % 计算近似春分点
sun_angle = 2*pi*(day_of_year-N0)/365.2422; % 计算太阳平均黄经
dec_angle = 0.3723+23.2567*sin(sun_angle)+0.1149*sin(2*sun_angle)-0.1712*sin(3*sun_angle)-0.758*cos(sun_angle)+0.3656*cos(2*sun_angle)+0.0201*cos(3*sun_angle);
dec_angle = deg2rad(dec_angle); % 计算赤纬角
time_angle = deg2rad((time-12)*15); % 计算太阳时角
hangle = asin(sin(latitude)*sin(dec_angle)+cos(latitude)*cos(dec_angle)*cos(time_angle));
hangle = rad2deg(hangle); % 计算太阳高度角
end

% 计算某一时刻被光伏板吸收的太阳辐射能量
function energy_3 = calculate_energy_3(year,month,day,time,incline_1,incline_2,k,energy_1)
% incline_1是光伏板水平仰角，incline_2是光伏板方位角
monthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
latitude = deg2rad(30.58333333333); % 转化纬度为弧度制
incline_1 = deg2rad(incline_1);
incline_2 = deg2rad(incline_2);
day_of_year = day;
for i = 1:month-1
    day_of_year = day_of_year + monthdays(i); % 计算年积日
end
N0 = 79.6764+0.2422*(year-1985)-floor((year-1985)/4); % 计算近似春分点
sun_angle = 2*pi*(day_of_year-N0)/365.2422; % 计算太阳平均黄经
dec_angle = 0.3723+23.2567*sin(sun_angle)+0.1149*sin(2*sun_angle)-0.1712*sin(3*sun_angle)-0.758*cos(sun_angle)+0.3656*cos(2*sun_angle)+0.0201*cos(3*sun_angle);
dec_angle = deg2rad(dec_angle); % 计算赤纬角
time_angle = deg2rad((time-12)*15); % 计算太阳时角
hangle = asin(sin(latitude)*sin(dec_angle)+cos(latitude)*cos(dec_angle)*cos(time_angle));
hangle = rad2deg(hangle); % 计算太阳高度角
cos_ang = cos((hangle+90)/180*pi);
energy_2 = energy_1 - k * energy_1*(6371*cos_ang+sqrt(6371^2*cos_ang.^2-6371^2+7371^2)); % 进入大气层到达地面时剩余能量(余弦定理求解)
cos_angle = sin(dec_angle)*(sin(latitude)*cos(incline_1)-cos(latitude)*sin(incline_1)*cos(incline_2))+cos(dec_angle)*(cos(latitude)*cos(incline_1)* ...
    cos(time_angle)+sin(latitude)*sin(incline_1)*cos(incline_2)*cos(time_angle)+sin(incline_1)*sin(incline_2)*sin(time_angle));
energy_3 = energy_2 .* cos_angle;
end