%%
% Experi�ncia 5: Cancelamento de Eco e RLS
% PSI3531 - Aplica��es de Filtragem Adaptativa
% Prof. 
% Matheus Bordin Gomes - 9838028

clear; close all; clc;

%%
% Dados do problema

% N�mero de amostras
N = 3500;
n = 1:N;
% Sinal x
var_x = 1;
% Sinal v
var_v = 0.01;
% Sistema
M = 100;
h = randn(M,1);
h = h/norm(h);

%%
% Execu��o do NLMS
mu = 0.5;
eps = 1e-3;
L = 100;
EMSE_nlms = zeros(N,1);
for i = 1:L
    x = sqrt(var_x)*randn(N,1);
    v = sqrt(var_v)*randn(N,1);
    [ e, ~, ~ ] = NLMS(x, v, h, mu, eps, h);
    EMSE_nlms = EMSE_nlms +(e.^2-var_v)/L;
end
fprintf('EMSE em regime do NLMS = %.5f\n', mean(EMSE_nlms(3000:end)));

%%
% Execu��o do RLS
lambda = 0.999;
EMSE_rls = zeros(N,1);
for i = 1:L
    x = sqrt(var_x)*randn(N,1);
    v = sqrt(var_v)*randn(N,1);
    [ e, ~, ~ ] = RLS(x, v, h, lambda, h);
    EMSE_rls = EMSE_rls +(e.^2-var_v)/L;
end
fprintf('EMSE em regime do RLS = %.5f\n', mean(EMSE_rls(3000:end)));

%%
% Plota o EMSE
figure();
semilogy(n, EMSE_nlms, n, EMSE_rls);
title('EMSE');
legend('NLMS', 'RLS');
xlabel('n (amostras)');
