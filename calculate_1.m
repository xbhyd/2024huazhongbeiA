% 从早6点到晚7点每隔半小时的时间序列
time = 6+(114.31666667-120)/15:0.5:19+(114.31666667-120)/15;
angle = calculate_height(2023,5,23,time,30.58333333333);
cos_angle = cos((angle+90)/180*pi);

% 计算穿越大气层的长度（余弦定理），经检验，选用算出的正值结果
s = 6371.*cos_angle+sqrt(6371^2*cos_angle.^2-6371^2+7371^2);

% 计算减少的太阳辐射
delta_x = 1334-[21.132 64.906 155.472 265.660 372.830 499.623 587.170 643.019 688.302 727.547 736.604 747.170 750.189 754.717 762.264...
    744.151 733.585 698.868 653.585 605.283 546.415 464.906 377.358 271.698 125.283 48.302 10.566];

% 计算衰减系数k
k = mean(delta_x/1334./s);
k_s = std(delta_x/1334./s); % 计算k的中误差

% 计算穿过了大气层但没有射入光伏板时剩余的太阳辐射
energy_1 = [1405 1394 1378 1353 1334 1316 1308 1315 1330 1350 1372 1392];
energy_3 = zeros(3,12); % 最大辐射结果
max_time = zeros(3,12); % 辐射结果最大对应的武汉地区的时间
for i=1:12
    for j=1:3
        f = @(t)calculate_energy_3(t,i,energy_1(1,i),j*20,k);
        options = optimset('Display','iter');
        [x_max, fval_max] = fminbnd(@(t) -f(t), 6, 18, options);
        energy_3(j,i) = -fval_max;
        max_time(j,i) = x_max;
    end
end
max_time_Beijing = max_time-(114.31666667-120)/15; % 辐射结果最大对应的北京地区的时间

% 接下来算每天1平方米受到的太阳直射辐射总能量（日出和日落时太阳高度角为0）
I = zeros(3,12);
for i=1:12
    for j=1:3
        f2 = @(t2)calculate_height(2025,i,15,t2,30.58333333333);
        % 计算日出日落时间
        t_min = fzero(f2,[3,8.5]);
        t_max = fzero(f2,[16,21]);
        f = @(t3)calculate_energy_3(t3,i,energy_1(1,i),j*20,k);
        I(j,i) = 3600 * integral(f, t_min, t_max); % 将单位转化为J每平方米，由于给定1平方米光伏板，则1天光伏板接受的能量数值上与此相等
    end
end

% 绘制图3
x_month = [1,2,3,4,5,6,7,8,9,10,11,12];
hold on;
plot(x_month,energy_3(1,:),"Color",'red','LineStyle','-');
plot(x_month,energy_3(2,:),"Color",'blue','LineStyle','-');
plot(x_month,energy_3(3,:),"Color",'green','LineStyle','-');
grid on;
xlabel('月份');
ylabel('最大太阳辐射');
legend('20°','40°','60°');

% 绘制图4
figure(2);
hold on;
plot(x_month,I(1,:),"Color",'red','LineStyle','-');
plot(x_month,I(2,:),"Color",'blue','LineStyle','-');
plot(x_month,I(3,:),"Color",'green','LineStyle','-');
grid on;
xlabel('月份');
ylabel('每月15日太阳辐射');
legend('20°','40°','60°');

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
function energy_3 = calculate_energy_3(t,month,energy_1,incline,k)
monthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
day_of_year = 15;
for i = 1:month-1
    day_of_year = day_of_year+monthdays(i); % 计算年积日
end
N0 = 79.6764+0.2422*(2025-1985)-floor((2025-1985)/4); % 计算近似春分点
sun_angle = 2*pi*(day_of_year-N0)/365.2422; % 计算太阳平均黄经
dec_angle = 0.3723+23.2567*sin(sun_angle)+0.1149*sin(2*sun_angle)-0.1712*sin(3*sun_angle)-0.758*cos(sun_angle)+0.3656*cos(2*sun_angle)+0.0201*cos(3*sun_angle);
dec_angle = deg2rad(dec_angle); % 计算赤纬角
time_angle = deg2rad((t-12)*15); % 计算太阳时角
cos_ang = cos((calculate_height(2025,month,15,t,30.58333333333)+90)/180*pi);
energy_2 = energy_1 - k * energy_1*(6371*cos_ang+sqrt(6371^2*cos_ang.^2-6371^2+7371^2)); % 进入大气层到达地面时剩余能量(余弦定理求解)
cos_angle = sin((30.58333333333-incline)/180*pi)*sin(dec_angle)+cos((30.58333333333-incline)/180*pi)*cos(dec_angle)*cos(time_angle); % 太阳光与光伏板法线夹角的余弦
energy_3 = energy_2 .* cos_angle; % 计算余弦损失剩下的能量
end