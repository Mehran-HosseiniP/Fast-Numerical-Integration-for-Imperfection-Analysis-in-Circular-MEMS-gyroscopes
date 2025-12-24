clear variables
close all

% Material properties (for silicon)
E = 170e+09; % Young's modulus (Pa)
rho = 2329; % mass density (kg/m^3)

% Ring dimensions (in m)
R = 1500e-06; % mean radius
w0 = 100e-06; % nominal width
h = 30e-06; % structural thickness

% Imprefection
dw = 200e-09; % peak-to-peak width variation
phi = pi/6; % max/min angle

syms theta real
w(theta) = 2*R - 2*sqrt((R-w0/2)^2 + dw^2/4 + dw*(R-w0/2)*cos(theta - phi));

figure,
set(gcf,'color','w')
fplot(w(theta*pi/180)*1e+06,[0, 360],'LineWidth',2)
hold on
fplot(w0*1e+06,[0, 360],'r--','LineWidth',2)
xticks(0:30:360)
grid on
set(gca,'LineWidth',2)
legend({'Varying width','Nominal width'},'Interpreter','latex','FontSize',13)
xlabel('$\theta(^\circ)$','Interpreter','latex','FontSize',14)
ylabel('$w(\theta)$ (um)','Interpreter','latex','FontSize',14)

% cross-section area
A = w*h;

% cross-section 2nd moment of area
I = w^3*h/12;

N = 8; % no. of intervals

theta_gl = GaussLegendre_discretization(0,2*pi,N);

%% Modal properties

% Mode number
n = 2;

% Mass matrix
A_gl = double(subs(A, theta, theta_gl));
int_A0 = legpoly_int(theta_gl, A_gl, 0); % zeroth harmonic
int_A1 = legpoly_int(theta_gl, A_gl, 2*n); % 2n-th harmonic

M = 0.5*rho*R*((n^2+1)/n^2)*real(int_A0)*eye(2) + ...
    0.5*rho*R*((n^2-1)/n^2)*[real(int_A1), imag(int_A1);
                             imag(int_A1), -real(int_A1)];

% Stiffness matrix
I_gl = double(subs(I, theta, theta_gl));
int_I0 = legpoly_int(theta_gl, I_gl, 0); % zeroth harmonic
int_I1 = legpoly_int(theta_gl, I_gl, 2*n); % 2n-th harmonic

K = 0.5*(E/R^3)*(n^2-1)^2*real(int_I0)*eye(2) + ...
    0.5*(E/R^3)*(n^2-1)^2*[real(int_I1), imag(int_I1);
                             imag(int_I1), -real(int_I1)];

%% ************************ Eigenvalue Problem ****************************

V_ref = [sin(2*n*phi), -sin(2*n*phi);
        (1-cos(2*n*phi)),1+cos(2*n*phi)];
V_ref = [V_ref(:,1)/norm(V_ref(:,1),2), V_ref(:,2)/norm(V_ref(:,2),2)];
% if det(V_ref)<0
%     V_ref = V_ref(:,[2, 1]);
% end
psi_ref = mode_rot(V_ref);
% V_ref = [cos(psi_ref), -sin(psi_ref);
%          sin(psi_ref), cos(psi_ref)];

[V,D] = eig(K, M);

% Natural frequencies
% freqs = sqrt([D(1,1),D(2,2)]);
% freq_split = abs(freqs(1) - freqs(2))/(2*pi);
% disp('Natural Frequencies (kHz):');
% display(freqs/(2000*pi));
% disp('Frequency Split (Hz):');
% display(freq_split);

% Mode shape rotation
V = [V(:,1)/norm(V(:,1),2), V(:,2)/norm(V(:,2),2)];
% if det(V)<0
%     V = V(:,[2, 1]);
% end
% pssi = mode_rot(V);
% V = [cos(pssi), -sin(pssi);
%          sin(pssi), cos(pssi)];

Vm{N} = (V*V')*(eye(2) - V_ref*V_ref');
Rho(N) = V_compare(V_ref, V);
psi_p{1} = [atan2(V(2,1),V(1,1)), atan2(-V(2,1),-V(1,1))];
psi_p{2} = [atan2(V(2,2),V(1,2)), atan2(-V(2,2),-V(1,2))];
psi_q{1} = [atan2(-V(1,1),V(2,1)), atan2(V(1,1),-V(2,1))];
psi_q{2} = [atan2(-V(1,2),V(2,2)), atan2(V(1,2),-V(2,2))];

[~,ind11] = min(abs(psi_p{1}));
psi11 = psi_p{1}(ind11);
[~,ind12] = min(abs(psi_p{2}));
psi12 = psi_p{2}(ind12);
[~,ind21] = min(abs(psi_q{1}));
psi21 = psi_q{1}(ind21);
[~,ind22] = min(abs(psi_q{2}));
psi22 = psi_q{2}(ind22);
if abs(psi11)<abs(psi21)
    psi1 = psi11;
    psi2 = psi22;
    % wn = [freqs(1), freqs(2)];
else
    psi1 = psi12;
    psi2 = psi21;
    % wn = [freqs(2), freqs(1)];
end

%%
psi1 = atan2d(V(2,1), V(1,1));
psi2 = atan2d(-V(2,1), -V(1,1));

disp('Mode shape rotation (deg) - mod pi/4:');
display([psi1/2, psi2/2]);


