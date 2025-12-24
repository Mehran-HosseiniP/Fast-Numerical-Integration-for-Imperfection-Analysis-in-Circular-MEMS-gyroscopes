clearvars
close all

syms x
f = (1+x+x^2+x^3);
w = 6;
F = f*cos(w*x);

% exact value
I_exact = double(int(F,x,-1,1));

% Gauss-Legendre integration
for i = 1:8
t_g = GaussLegendre_discretization(-1,1,i);
x_g = double(subs(f,x,t_g));
I_g = real(legpoly_int(t_g,x_g,w));
e_g(i) = abs((I_g-I_exact));
nodes_g(i) = length(t_g);
end
% fprintf('Gauss-Legendre error:%.4f   no. of nodes:%0.4f\n',e_g,numel(t_g));

% Simpson 1/3 integration
for i = 1:8
t_s = linspace(-1,1,nodes_g(i)+1);
x_s = double(subs(F,x,t_s));
I_s = simpson13_int(t_s,x_s);
e_s(i) = abs((I_s-I_exact));
end
% fprintf('Simpson 1/3 error:%.4f   no. of nodes:%0.4f\n',e_s,numel(t_s));

for i = 1:8
t_p = linspace(-1,1,nodes_g(i));
x_p = double(subs(F,x,t_p));
I_p = trapezoidal_int(t_p,x_p);
e_p(i) = abs((I_p-I_exact));
end

figure(1)
set(gcf,'color','w')
plot(nodes_g, log10(e_g),'o-','color','b','MarkerSize',13,'linewidth',1.5)
hold on
plot(nodes_g+1, log10(e_s),'s--','color','r','MarkerSize',13,'linewidth',1.5)
hold on
plot(nodes_g, log10(e_p),'*--','color','g','MarkerSize',13,'linewidth',1.5)
xlabel('\#Nodes','interpreter','latex','fontsize',14)
ylabel('$\log_{10}$ of Error','interpreter','latex','fontsize',14)
legend({'Gauss-Legendre','Simpson 1/3','Trapezoidal'},'interpreter','latex','fontsize',14)