clear;
close all;

xin = csvread('x0.dat');
% x_s = reshape(xin', [], 2);
% csvwrite("E:\DATA\xin1.dat",x_s(:,1));
% csvwrite("E:\DATA\xin2.dat",x_s(:,2));


x_l = 600;         % �źŽ�ȡ����
b_l = 100;          % ���˲��㳤��

FS = 200e3;     % 200kHz����

N = 8;          % ��Ԫ����
R = 0.24;       % Բ��뾶
BETA = 45*pi/180;   % ��Ԫ�Ƕȼ��

C = 1500;       % ����

theta = (-180:30:180)*pi/180; 
x_d = zeros(N,x_l+b_l*2);
delay = zeros(N,length(theta));
D = zeros(length(theta),1);
for t = 1:length(theta)-1
    for i = 1:N
        delay(i,t) = FS*R*cos(theta(t)-(i-1)*BETA)/C;
        x_d(i,:) = [zeros(1,b_l+round(delay(i,t))), xin(i,1:x_l), zeros(1,b_l-round(delay(i,t)))];
    end
    w = [1 1 1 1 1 1 1 1];      % ���ȼ�Ȩ
    w = CRS(w,t-1);
    x_dd = diag(w)*x_d;
    ccs_d = reshape(xx_d,[],8)';
    y = sum(x_dd);
    D(t) = sum(abs(y));    
end
[~ , n] = max(D);
v1 = D(n-1);
v2 = D(n);
v3 = D(n+1);

u = -180+(n-1)*30;                         % �������������������������
A_max = v1-v3;
B_max = v1+v3-2*v2;
bizhi = A_max/B_max;
x = u+bizhi/2*30;                           % x = u+bizhi/2*�ನ��֮��Ŀ��;
fprintf('����ֵ��%f \n',x);