function x = GaussLegendre_discretization(lower_bound, upper_bound, N)

% Generating nodes for Gauss-Legendre quadrature over [lower_bound,
%upper_bound] using N equal intervals

% Zeros of Legendre polynomial P_4(x)
zer1 = sqrt(3/7-(2/7)*sqrt(6/5));
zer2 = sqrt(3/7+(2/7)*sqrt(6/5));

x_dummy = linspace(lower_bound, upper_bound, N+1);
A = [1 + zer2, 1 - zer2;
     1 + zer1, 1 - zer1;
     1 - zer1, 1 + zer1;
     1 - zer2, 1 + zer2]/2;
x = zeros(4*N, 1);

for i = 1:N
    x(4*i-3:4*i, 1) = A*[x_dummy(i);x_dummy(i+1)];
end

end