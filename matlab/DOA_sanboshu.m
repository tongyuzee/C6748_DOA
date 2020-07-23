% ====================== 方位估计实验 三波束比幅法 ======================= %
% ============================== 8元圆阵 ================================ %
clear ;
close all;
% 导入数据
load('data_831.mat');
x0(1,:) = Data1_AI_0___U';
x0(2,:) = Data1_AI_1___U';
x0(3,:) = Data1_AI_2___U';
x0(4,:) = Data1_AI_3___U';
x0(5,:) = Data1_AI_4___U';
x0(6,:) = Data1_AI_5___U';
x0(7,:) = Data1_AI_6___U';
x0(8,:) = Data1_AI_7___U';
csvwrite("x0.dat",x0(:,146000:147000-1));
% 基阵信息
N = 8;
R = 0.24;
BETA = 45*pi/180;  % 阵元角度间隔
% 信号信息
% F = 5700;        % 5.7khz信号
FS = 200e3;     % 200kHz采样
% T = 100/F;     % 100个CYCLE
% T0 = 1;        % 1s周期
% NS = T*FS;     % 脉宽内的采样点
C = 1500;        % 声速
THETA_D = 30;    % 多波束的波束指向角增量,必须是360的约数
THETA_T = (-180:THETA_D:180)*pi/180;                    % 假定信号来得方向
x = [zeros(N,300) x0(:,146000:147000) zeros(N,300)]; % 信号长度（根据实际情况调整）
M = length(x(1,:)); 
figure(1)
for i=1:N
    subplot(2,4,i) 
    plot(x(i,:)); 
    axis([0 M-1 -5 5])
    title(strcat(num2str(i),'号阵元'));
end
%% 多波束形成/各波束幅度响应
for t = 1:length(THETA_T)-1                    % 形成多波束，波束号t
    % 在假定入射方向THETA_T(t)时，各阵元补时延差，使其同相叠加，提前为正
    for i = 0:N-1                              % 阵元号i
        delay_T(i+1) = FS*R*cos(THETA_T(t)-i*BETA)/C;                      % 相对相位中心的延迟为-，提前为+
        x_d(i+1,:) = interp1(0:M-1,x(i+1,:),(0:M-1)-delay_T(i+1),'cubic'); % 提前量推迟回去 
        x_dd(i+1,:) = [zeros(1,300+round(delay_T(i+1))), x0(i+1,146000:147000), zeros(1,300-round(delay_T(i+1)))];
%         figure(10);plot(x_d(i+1,:));hold on; plot(x(i+1,:))
        
        % vq=interp1（x，v，xq）使用线性插值返回一维函数在特定查询点的插值值，默认方法是“linear”
        % x:采样点0:M-1，v:相应的值x(i+1,:)，xq:查询点的坐标(0:M-1)-delay_T(i+1)
        % vq = interp1(x,v,xq,method)指定了插值方法，cubic：保形分段三次插值,查询点上的插值值是基于相邻网格点上的值的保形分段三次插值
    end
%     w = [1 1 1 0 0 0 1 1];      % 0°指向角时的加权向量
    w = [1 1 1 1 1 1 1 1];      % 均匀加权
    % 调用循环移位子函数
    w = CRS(w,t-1);
    x_d = diag(w)*x_d;
    y = sum(x_d);
    D(t) = sum(abs(y));         % 在假定入射方向THETA_T(t)时，该波束输出幅度响应
end
[Dmax1,t1] = maxi(D,1:length(THETA_T)-1);   % 寻找最大值和对应下标
figure(2) 
stem(1:length(THETA_T)-1,D,'*'); 
title('波束输出幅度响应');
xlabel('波束号t');
ylabel('各波束幅度响应');
fprintf('最大输出波束：%d \n',t1);
% 三个波束的输出幅值
v1=D(t1-1);                                 % 左侧相邻输出
fprintf('左侧波束幅值：%d .\n',v1);
v2=D(t1);                                   % 中心波束最大输出
fprintf('中心波束幅值：%d .\n',v2);
v3=D(t1+1);                                 % 右侧相邻输出
fprintf('右侧波束幅值：%d .\n',v3);
%% 三波束比幅法
u = -180+(t1-1)*30;                         % 最大输出波束的主极大方向的坐标
A_max = v1-v3;
B_max = v1+v3-2*v2;
bizhi = A_max/B_max;
x = u+bizhi/2*30;                           % x = u+bizhi/2*多波束之间的宽度;
fprintf('测向值：%f \n',x);


