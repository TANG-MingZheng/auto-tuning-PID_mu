function [Kp,Ti,Td] = Function_ZNPID_mu(mu,Kc,Tc)
%Function_ZNPID_mu ��������ʽ��C(s) = Kp(1 + 1/(Ti*s) + Td*s^mu)
%   ����������ٽ��񵴲����ͷ����׵Ľ״Σ����ط����� PID�� ����������������
%   Kc  �ٽ�������
%   Tc  �ٽ�������
p = (0.6*Kc);
q = (0.6*Kc*(pi/4-1/pi));
m = q / p;
C = cos(mu*pi/2);
S = sin(mu*pi/2);
a = ((m*C-S)*4);
b = (4*m);
c = 1;
syms x;
equ = a*x^(1+mu) + b*x + c;
% dif = 4*(mu+1)*(m*C-S)*x^mu + 4*m;
% figure;
% hold on;
% fplot(equ);
% fplot(dif);
% grid on;
% hold off;
xv = double(vpasolve(equ, x));
if length(xv) == 1
    X = xv;
elseif length(xv) > 1
    X = max(xv);
else
    error('No soluation !');
end
alpha = 0.6/(1+C*X^mu);
beta  = 2*X/pi;
gamma = (X/(2*pi))^mu;
Kp = alpha*Kc;
Ti = beta*Tc;
Td = gamma*Tc^mu;
end

