function [m,i] = maxi(x,n)
%Ѱ�����ֵ������±�
%x������У�n�������е��±�
%max�������ֵ��i�������ֵ�±�
m = max(x);
for i = n(1):length(n)+n(1)-1
    if(x(i)==m)
        break;
    end
end
