clear
clc
close all

mu_0 = pi*4e-7;
eps_0 = 8.854e-12;
eps_r1 = 2; %I choose
eps_r2 = 3;
c = 1/sqrt(mu_0*eps_0);

f = 20e6; %I choose
lambda_0 = c/f;
lambda_1 = lambda_0/sqrt(eps_r1);
lambda_2 = lambda_0/sqrt(eps_r2);
lambda_3 = lambda_0/sqrt(eps_r1);

deltax = 0.02; %I choose
xmax = 175; %I choose
x = 0:deltax:xmax;

tmax = xmax/c; %I choose
S = 1; %I choose
deltat = S*deltax/c;
t = 0;
d1 = lambda_1/4;
d2 = lambda_2/2;
d3 = lambda_3/4;

eps = ones(length(x), 1) .* eps_0;
s1 = round(((xmax/2) - ((d1+d2+d3)/2))/0.02); %I choose
s2 = s1 + round(d1/deltax); %makes dielectric slab of thickness d
s3 = s2 + round(d2/deltax);
s4 = s3 + round(d3/deltax);
eps(s1:s2) = eps_0*eps_r1;
eps(s2:s3) = eps_0*eps_r2;
eps(s3:s4) = eps_0*eps_r1;


E = zeros(length(x),3);

x0 = 1;
x1 = 5; %I choose

t1 = 150/c; %I choose
n = 0;
Ncyc = 3; %number of oscillations pulsed out
E0max = 0;
E1max = 0;

while t<tmax
    E(2:end-1,3) = ((1./(mu_0*eps(2:end-1))) .* (deltat/deltax)^2 .* (E(3:end,2) - 2*E(2:end-1,2) + E(1:end-2,2))) + 2*E(2:end-1,2) - E(2:end-1,1);
    
    if t<Ncyc/f
        E(x0,3) = sin(2*pi*f*t);
    else
        E(x0,3) = 0;
    end

    n = n + 1;

    if mod(n,100) == 0
        plot(x, E(:,3)); hold on
        xline(x(s1), 'r--'); xline(x(s2),'r--'); xline(x(s3),'r--'); xline(x(s4),'r--'); hold off
        ylim([-1.5 1.5]); xlabel('x (m)'); ylabel('E')
        title(sprintf('f = %g MHz, t = %.0f ns', f/1e6, t*1e9))
        drawnow
    end

    E(:,1) = E(:,2);
    E(:,2) = E(:,3);
    t = t + deltat;
end