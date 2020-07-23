function [y]=CRS(x,n0)
% CRS：Circle Rigth Shift 循环右移
% n0为所要移动的位数
for i = 0:length(x)-1
    l = mod((i+n0),length(x));
    % b = mod(a,m)返回a除以m后的余数，其中a是被除数，m是除数。
    y(l+1)=x(i+1);
end
end