%%
% Experiência 5: Cancelamento de Eco e RLS
% PSI3531 - Aplicações de Filtragem Adaptativa
% Prof. 
% Matheus Bordin Gomes - 9838028

clear; close all; clc;

%%
% Dados do problema

% Sinal x
[x,Fsx] = audioread('conversa2.wav');
% Sinal v
[v,Fsv] = audioread('conversa1.wav');
% Sistema
load('respimp.mat')
h = ri256;
M = length(h);
% Tamanho dos vetores
N = length(x);
n = 1:N;

%%
% Execução do NLMS
mu = 0.01;
eps = 0.1;
tic
[ e_nlms, delta_w_nlms, y_nlms, runtime ] = NLMS(x, v, h, mu, eps, h);
nlms_avg_time = runtime;
fprintf('NLMS avarage time = %.5f\n', nlms_avg_time);
ERLE_nlms = 10*log10(y_nlms.^2./e_nlms.^2);
MSD_nlms = sum(delta_w_nlms.^2,2);

%%
% Execução do RLS
lambda = 0.999999;
[ e_rls, delta_w_rls, y_rls, runtime ] = RLS(x, v, h, lambda, h);
rls_avg_time = runtime;
fprintf('RLS avarage time = %.5f\n', rls_avg_time);
ERLE_rls = 10*log10(y_rls.^2./e_rls.^2);
MSD_rls = sum(delta_w_rls.^2,2);

%%
% Comparação entre e[n] e v[n]
figure();
subplot(3,1,1);
plot(n, v);
title('v[n]');
xlabel('n (amostras)');
subplot(3,1,2);
plot(n, e_nlms);
title('e[n] (NLMS)');
xlabel('n (amostras)');
subplot(3,1,3);
plot(n, e_rls);
title('e[n] (RLS)');
xlabel('n (amostras)');

%%
% Plota o ERLE
figure();
subplot(2,1,1);
plot(n, ERLE_nlms);
title('ERLE (NLMS)');
xlabel('n (amostras)');
subplot(2,1,2);
plot(n, ERLE_rls);
title('ERLE (RLS)');
xlabel('n (amostras)');

%%
% Plota o MSD
figure();
semilogy(n, MSD_nlms, n, MSD_rls);
title('MSD');
legend('NLMS', 'RLS');
xlabel('n (amostras)');
