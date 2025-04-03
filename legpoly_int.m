function [y] = legpoly_int(t,x,w)
% Numerical integration of x(t)*exp(1j*w*t) using Legendre polynomials
% Last modified: Oct 2, 2021

%*********************************************************************
% Number of nodes should be 4*N (4,8,...)
N = floor(numel(t)/4);

%*********************************************************************
% Zeros of the Legendre polynomial P4(x)
zer1 = sqrt(3/7-(2/7)*sqrt(6/5));
zer2 = sqrt(3/7+(2/7)*sqrt(6/5));

%*********************************************************************
% Gauss-Legendre weights
weigh1 = (18+sqrt(30))/36;
weigh2 = (18-sqrt(30))/36;

%*********************************************************************
% Nodal values of Legendre polynomials
A = [1.0000   -0.8611    0.6123   -0.3047
    1.0000   -0.3400   -0.3266    0.4117
    1.0000    0.3400   -0.3266   -0.4117
    1.0000    0.8611    0.6123    0.3047];

%**********************************************************************
% Integration
y = 0;

for i = 1:N
    da = (t(4*i-1)-t(4*i))/(zer1-zer2);
    a_avg = (zer1*t(4*i)-zer2*t(4*i-1))/(zer1-zer2);
    w1 = w*da;
    
    % Integrals of legendre polynomial P0(t),...,P3(t)
    
    if w==0
       Jp = [2 0 0 0];
    else
       for kk = 1:4
           Jp(kk) = sqrt(2*pi/w1)*(1j)^(kk-1)*besselj(kk-0.5,w1);
       end
    end
    
    % Coefficients (approximated by Gauss–Legendre quadrature)
    
    xvec = [weigh2*x(4*i-3);weigh1*x(4*i-2);weigh1*x(4*i-1);weigh2*x(4*i)];
    alpha = zeros(4,1);
    
    for kj = 1:4
        alpha(kj) = ((2*kj-1)/2)*dot(xvec,A(:,kj));
    end
    
    y = y + da*exp(1j*w*a_avg)*dot(alpha,Jp);
    
end

end
