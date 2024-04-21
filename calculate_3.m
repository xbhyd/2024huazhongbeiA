% 模型设置
% 适应度函数的函数句柄
fitnessfcn = @fun;
% 变量个数
nvars = 2;
% 约束条件形式1：下限与上限（若无取空数组[]）
% lb<= X <= ub
lb = [0, -90];
ub = [90, 90];
A = [];
b = [];
% 约束条件形式3：线性等式约束（若无取空数组[]）
% Aeq*X == beq
Aeq = [];
beq = [];
% 求解器设置
% 最优个体系数paretoFraction
% 种群大小populationsize
% 最大进化代数generations
% 停止代数stallGenLimit
% 适应度函数偏差TolFun
% 函数gaplotpareto：绘制Pareto前沿 
options = gaoptimset('paretoFraction', 0.3, 'populationsize', 200, 'generations', 300, 'stallGenLimit', 200, 'TolFun', 1e-10, 'PlotFcns', @gaplotpareto);
% 主求解
[x, fval] = gamultiobj(fitnessfcn, nvars, A, b, Aeq, beq, lb, ub, options);

% 结果提取
% 因为gamultiobj是以目标函数分量取极小值为目标，
% 因此在y=Fun(x)里取相反数的目标函数再取相反数画出原始情况
plot(-fval(:, 1)/365, -fval(:, 2)/365, 'pr');
xlabel('f_1(x)');
ylabel('f_2(x)');
title('Pareto front');
grid on;

% 绘制仰角与方位角对应关系图
incline_1 = x(:,1);
incline_2 = x(:,2);
figure(2);
plot(incline_2, incline_1, 'pr');
xlabel('方位角');
ylabel('水平仰角');
grid on;

function y = fun(x)
    incline_1 = x(1);
    incline_2 = x(2);
    y(1) = -calculate_year_energy(incline_1, incline_2);
    y(2) = -calculate_year_time(incline_1, incline_2);
end

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
            f2 = @(t2)calculate_energy_3(2025, i, j, t2, incline_1, incline_2, k, energy_1(i),0);
            energy_of_year = energy_of_year + 3600 * integral(f2, t_min, t_max);
        end
    end
end

% 计算一年内符合最理想情况的时间长度
function time_of_year = calculate_year_time(incline_1, incline_2)
    monthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
    energy_1 = [1405 1394 1378 1353 1334 1316 1308 1315 1330 1350 1372 1392];
    time_of_year = 0;
    k = 0.000401466772612056;
    for i=1:12
        for j=1:monthdays(i)
            disp(num2str(i)+" " +num2str(j))
            % 计算太阳辐射上午到达150W/m2的时刻
            f3 = @(t3)calculate_energy_3(2025, i, j, t3, incline_1, incline_2, k, energy_1(i), 150);
            t_min = fsolve(f3,7);
            % 计算太阳辐射下午到达100W/m2的时刻
            f4 = @(t4)calculate_energy_3(2025, i, j, t4, incline_1, incline_2, k, energy_1(i), 100);
            t_max = fsolve(f4,18);
            time_of_year = time_of_year + (t_max - t_min);
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
function energy_3 = calculate_energy_3(year,month,day,time,incline_1,incline_2,k,energy_1,t)
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
energy_3 = energy_2 .* cos_angle - t;
end