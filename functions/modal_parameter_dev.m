function [omega_dev, zeta_dev] = modal_parameter_dev(Lambda1, Lambda2)

omega1 = abs(Lambda1);
omega2 = abs(Lambda2);

zeta1 = real(Lambda1) ./ omega1;
zeta2 = real(Lambda2) ./ omega2;

omega_dev = dev(omega1, omega2);
zeta_dev = dev(zeta1, zeta2);

end