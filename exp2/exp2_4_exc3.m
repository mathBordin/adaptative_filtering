%%
% Experiência 3: Identificação de Sistemas
% PSI3531 - Prof. Vítor
% Matheus Bordin Gomes - 9838028
clear;clc;close all;
%#ok<*NOPTS>

%%
% Criação dos parâmetros

N = 50000;
n = 1:N;
M = 5;
% Sinal x[n]
var_x = 1;
a=0.9;
Hb_in=(1-a^2);
Ha_in=[1 -a];
% Sinal v[n]
var_v = 1e-2;
% Sistema
H = [1 -0.5 0 0 0];
% Coeficientes do filtro adaptativo
W = zeros(N,M);
% Passo 
mu = 0.005;
% W ótimo
R_phi = zeros(M);
for i=1:M
    for j=1:M
        R_phi(i,j)=(a^abs(i-j))*(1-a^2);
    end
end
p = [(1-a^2)*H(1); (1-a^2)*H(2); 0; 0; 0];
Wo = (R_phi\p);


%%
% Execuções do algoritmo

K = 100;
% Diferença entre os coeficientes W em relação aos de H
delta_W = zeros(N,M); 
% Erro médio quadrático
MSE = zeros(N,1);
% Erro médio quadrático em excesso
EMSE = zeros(N,1);
% Desvio médio quadrático
MSD = zeros(N,1);
% Erro
erro = zeros(N,1);
% Execuções
for i = 1:K
    % Sinal x[n]
    x = sqrt(var_x)*randn(N,1);
    x = filter(Hb_in,Ha_in,x);
    % Sinal v[n]
    v = sqrt(var_v)*randn(N,1);
    % LMS
    [ W, e, y_est, emse, msd] = LMS_under_test(x, v, M, N, H, mu, H);
    % Parâmetros para avaliar o aprendizado do filtro
    MSE = MSE+e.^2;
    EMSE = EMSE + emse;
    MSD = MSD + msd;
    erro = erro + e;
end
MSE=MSE/K;
EMSE=EMSE/K;
MSD=MSD/K;
erro = erro/K;

% Cálculo das curvas teóricas
K_teo = H'*H;
sigma_0 = var(erro);

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

% Plota parâmetros de aprendizado
figure(1);
plot(n,MSE,n,MSE_teo);
title('Erro quadrático médio - MSE');
xlabel('Iterações');
ylabel('MSE');
legend('experimental','teórico');

figure(2);
plot(n,EMSE,n,EMSE_teo);
title('Erro quadrático médio em excesso - EMSE');
xlabel('Iterações');
ylabel('EMSE');
legend('experimental','teórico');

figure(3);
plot(n,MSD,n,MSD_teo);
title('Desvio médio quadrático - MSD');
xlabel('Iterações');
ylabel('MSD');
legend('experimental','teórico');

figure()
plot(1:N,W(:,1),1:N,W(:,2),1:N,W(:,3),1:N,W(:,4),1:N,W(:,5));
title('Evolução dos coeficientes do filtro adaptativo');
xlabel('Iterações');
ylabel('Evolução dos coeficientes');
legend('w1','w2', 'w3', 'w4', 'w5');