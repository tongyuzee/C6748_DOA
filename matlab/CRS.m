function [y]=CRS(x,n0)
% CRS��Circle Rigth Shift ѭ������
% n0Ϊ��Ҫ�ƶ���λ��
for i = 0:length(x)-1
    l = mod((i+n0),length(x));
    % b = mod(a,m)����a����m�������������a�Ǳ�������m�ǳ�����
    y(l+1)=x(i+1);
end
end