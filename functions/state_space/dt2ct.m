function [A_c, B_c] = dt2ct(A_d, B_d, dt)
% [A_c, B_c] = DT2CT(A_d, B_d, dt)
% Computes continuous-time state matrix A_c and continuous time input
% matrix B_c.
%
% Arguments:
% A_d: discrete-time state matrix
% B_d: discrete-time input matrix
% dt: time increment

A_c = logm(A_d)/dt;
A_c = round(A_c,12);

B_c = (A_c\A_d-inv(A_c))\B_d;
B_c = round(B_c, 14);

