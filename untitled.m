clear
clc
close all

mu_0 = pi*4e-7;
eps_0 = 8.854e-12;
eps_r = ; %I choose
c = 1/sqrt(mu_0*eps_0);

f = ; %I choose
lambda = c/f;

deltax = ; %I choose
xmax = ; %I choose
x = 0:deltax:xmax;

tmax = ; %I choose
S = ; %I choose
deltat = S*deltax/c;
t = 0;

eps = ones(length(x), 1) .* eps_0;
s1 = ; %I choose
eps(s1:end) = eps(s1:end)*eps_r;

E = zeros(length(x), 3);

x0 = 1;
x1 = ; %I choose

t1 = ; %I choose

while t<tmax
    E(2:end-1,3) = ;
end