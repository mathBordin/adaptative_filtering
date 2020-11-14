%%
% Experiência 5: Cancelamento de Eco e RLS
% PSI3531 - Aplicações de Filtragem Adaptativa
% Prof. 
% Matheus Bordin Gomes - 9838028

clear; close all; clc;

%%
% Dados do problema

% Número de amostras
N = 5000;
n = 1:N;
% Sinal x
var_x = 1;
a=0.9;
Hb_in=sqrt(1-a^2);
Ha_in=[1 -a];
% Sinal v
var_v = 0.01;
% Sistema
M = 100;
h = randn(M,1);
h = h/norm(h);

%%
% Execução do NLMS
mu = 0.6;
eps = 1e-4;
L = 100;
EMSE_nlms = zeros(N,1);
nlms_avg_time = 0;
for i = 1:L
    x = sqrt(var_x)*randn(N,1);
    x = filter(Hb_in,Ha_in,x);
    v = sqrt(var_v)*randn(N,1);
    [ e, ~, ~, runtime ] = NLMS(x, v, h, mu, eps, h);
    nlms_avg_time = nlms_avg_time + runtime/L;
    EMSE_nlms = EMSE_nlms +(e.^2-var_v)/L;
end
fprintf('EMSE em regime do NLMS = %.5f\n', mean(EMSE_nlms(4500:end)));
fprintf('NLMS avarage time = %.5f\n', nlms_avg_time);

%%
% Execução do RLS
lambda = 0.99;
EMSE_rls = zeros(N,1);
rls_avg_time = 0;
for i = 1:L
    x = sqrt(var_x)*randn(N,1);
    x = filter(Hb_in,Ha_in,x);
    v = sqrt(var_v)*randn(N,1);
    [ e, ~, ~, runtime ] = RLS(x, v, h, lambda, h);
    rls_avg_time = rls_avg_time + runtime/L;
    EMSE_rls = EMSE_rls +(e.^2-var_v)/L;
end
fprintf('EMSE em regime do RLMS = %.5f\n', mean(EMSE_rls(3000:end)));
fprintf('RLS avarage time = %.5f\n', rls_avg_time);

%%
% Plota o EMSE
figure();
semilogy(n, EMSE_nlms, n, EMSE_rls);
title('EMSE');
legend('NLMS', 'RLS');
xlabel('n (amostras)');
