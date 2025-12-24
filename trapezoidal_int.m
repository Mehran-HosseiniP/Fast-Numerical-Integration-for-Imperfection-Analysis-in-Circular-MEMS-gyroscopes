function [y1] = trapezoidal_int(t,x)
%Trapezoidal integration
y1 = 0;
for n = 1:numel(t)-1
    y1 = y1 + (t(n+1)-t(n))*(x(n+1)+x(n))/2;
end
end