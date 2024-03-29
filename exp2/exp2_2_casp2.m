%%
% Experi�ncia 3: Identifica��o de Sistemas
% PSI3531 - Prof. V�tor
% Matheus Bordin Gomes - 9838028
clear;clc;close all;
%#ok<*NOPTS>

%%
% Cria��o dos par�metros

N = 2000;
n = 1:N;
M = 3;
% Sinal x[n]
var_x = 1;
% Sinal v[n]
var_v = 1e-2;
% Sistema
H = [1 -0.5 0.2];
% Coeficientes do filtro adaptativo
W = zeros(N,M);
% Passo 
mu = 0.005;

%%
% Execu��es do algoritmo

K = 100;
% Diferen�a entre os coeficientes W em rela��o aos de H
delta_W = zeros(N,M); 
% Erro m�dio quadr�tico
MSE = zeros(N,1);
% Erro m�dio quadr�tico em excesso
EMSE = zeros(N,1);
% Desvio m�dio quadr�tico
MSD = zeros(N,1);
% Erro
erro = zeros(N,1);
% Execu��es
for i = 1:K
    % Sinal x[n]
    x = sqrt(var_x)*randn(N,1);
    % Sinal v[n]
    v = sqrt(var_v)*randn(N,1);
    % LMS
    [ W, e, y_est, emse, msd] = LMS_under_test(x, v, M, N, H, mu, H);
    % Par�metros para avaliar o aprendizado do filtro
    MSE = MSE+e.^2;
    EMSE = EMSE + emse;
    MSD = MSD + msd;
    erro = erro + e;
end
MSE=MSE/K;
EMSE=EMSE/K;
MSD=MSD/K;
erro = erro/K;

% C�lculo das curvas te�ricas
K_teo = H'*H;
sigma_0 = var(erro);
R_phi = [var_x 0 0; 0 var_x 0; 0 0 var_x]; 

MSD_teo = zeros(N,1);
EMSE_teo = zeros(N,1);
MSE_teo = zeros(N,1);
    
MSD_teo(1) = trace(K_teo);
EMSE_teo(1) = trace(R_phi*K_teo);
MSE_teo(1) = EMSE_teo(1) + sigma_0;
for i = 2:N
    K_teo_next = K_teo-(mu*R_phi*K_teo)-(mu*K_teo*R_phi)+((mu^2)*(sigma_0^2)*R_phi);
    K_teo = K_teo_next;
    MSD_teo(i) = trace(K_teo);
    EMSE_teo(i) = trace(R_phi*K_teo);
    MSE_teo(i) = EMSE_teo(i) + sigma_0;
end

% Plota par�metros de aprendizado
figure(1);
plot(n,MSE,n,MSE_teo);
title('Erro quadr�tico m�dio - MSE');
xlabel('Itera��es');
ylabel('MSE');
legend('experimental','te�rico');

figure(2);
plot(n,EMSE,n,EMSE_teo);
title('Erro quadr�tico m�dio em excesso - EMSE');
xlabel('Itera��es');
ylabel('EMSE');
legend('experimental','te�rico');

figure(3);
plot(n,MSD,n,MSD_teo);
title('Desvio m�dio quadr�tico - MSD');
xlabel('Itera��es');
ylabel('MSD');
legend('experimental','te�rico');