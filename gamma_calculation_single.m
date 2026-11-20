clear
clc
close all

mu_0 = pi*4e-7;
eps_0 = 8.854e-12;
eps_r1 = 2; %I choose
c = 1/sqrt(mu_0*eps_0);

%f = 20e6; %I choose
lambda_0 = c/20e6;
lambda_1 = lambda_0/sqrt(eps_r1);

deltax = 0.02; %I choose
xmax = 175; %I choose
x = 0:deltax:xmax;

%tmax = xmax/c; %I choose
S = 1; %I choose
deltat = S*deltax/c;
t = 0;
d1 = lambda_1/2;

eps = ones(length(x), 1) .* eps_0;
s1 = round(((xmax/2) - ((d1)/2))/0.02); %I choose
s2 = s1 + round(d1/deltax); %makes dielectric slab of thickness d
eps(s1:s2-1) = eps_0*eps_r1;

E = zeros(length(x),3);

x0 = 1;
x1 = s2 + round(5/deltax); %I choose
tmax = (x(x1) + 2*min(x(s1), xmax - x(x1)))/c;

fsweep = [10e6 12e6 14e6 16e6 18e6 20e6 22e6 24e6 26e6 28e6 30e6];
gamma_fdtd = zeros(size(fsweep));

for k = 1:numel(fsweep)
    f = fsweep(k);
    E = zeros(length(x),3);
    t = 0;
    E0max = 0;
    E1max = 0;

    while t<tmax
        E(2:end-1,3) = ((1./(mu_0*eps(2:end-1))) .* (deltat/deltax)^2 .* (E(3:end,2) - 2*E(2:end-1,2) + E(1:end-2,2))) + 2*E(2:end-1,2) - E(2:end-1,1);
        E(x0,3) = sin(2*pi*f*t);
    
    if t>tmax -2/f
        E0max = max(E0max, abs (E(x0,3)));
        E1max = max(E1max, abs (E(x1,3)));
    end

    E(:,1) = E(:,2);
    E(:,2) = E(:,3);
    t = t + deltat;
    end

tau = E1max/E0max;
gamma_fdtd(k) = sqrt(1-tau^2);

end

eta0 = sqrt(mu_0/eps_0);
fa = linspace(5e6, 35e6, 500);
er = [eps_r1];
dd = [d1];
Z = eta0*ones(size(fa));
for m = 1:-1:1
    eta = eta0/sqrt(er(m));
    tb = tan(2*pi*fa*sqrt(er(m))*dd(m)/c);
    Z = eta*(Z + 1j*eta*tb)./(eta + 1j*Z.*tb);
end
gamma_th = abs((Z-eta0)./(Z+eta0));

figure
plot(fa/1e6, gamma_th, 'LineWidth', 1.5); hold on
plot(fsweep/1e6, gamma_fdtd, 'o', 'MarkerSize', 7, 'LineWidth', 1.5)
xlabel('Frequency (MHz)'); 
ylabel('|\Gamma|'); 
legend('Analytic', 'FDTD'); grid on