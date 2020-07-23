function [m,i] = maxi(x,n)
%寻找最大值和最大下标
%x输出序列，n输入序列的下标
%max返回最大值，i返回最大值下标
m = max(x);
for i = n(1):length(n)+n(1)-1
    if(x(i)==m)
        break;
    end
end
