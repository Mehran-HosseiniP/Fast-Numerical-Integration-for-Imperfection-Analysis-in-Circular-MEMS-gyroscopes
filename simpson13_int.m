function y = simpson13_int(t,x)
%Simpson's 1/3 integration
N = length(t);
y = 0;
i = 1;

while i+2<=N
    y = y + (t(i+2)-t(i))*(x(i+2)+4*x(i+1)+x(i))/6;
    i = i+2;
end

if i==N-1
   x_mid = (x(N)+x(N-1))/2;
   y = y + (t(N)-t(N-1))*(x(N)+4*x_mid+x(N-1))/6;
end

end