function [omega_dev, zeta_dev, Lambdas] = comp_models(A1, A2)

Lambda1 = eig(A1);
[~, ind] = sortrows(abs(Lambda1));
Lambda1 = Lambda1(ind);

Lambda2 = eig(A2);
[~, ind] = sortrows(abs(Lambda2));
Lambda2 = Lambda2(ind);

omegas1 = abs(Lambda1);
omegas2 = abs(Lambda2);

zetas1 = real(Lambda1) ./ omegas1;
zetas2 = real(Lambda2) ./ omegas2;

omega_dev = dev(omegas1, omegas2);
zeta_dev = dev(zetas1, zetas2);

Lambdas = [Lambda1, Lambda2];