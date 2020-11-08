%%
% Experiência 4: Rastreamento
% PSI3531 - Prof. Vítor
% Matheus Bordin Gomes - 9838028
clear;clc;close all;
%#ok<*NOPTS>

%%
% Exercício 1)

%%
% Parâmetros iniciais
N = 2000;
n = 1:N;
M = 2;
% Sinal x[n]
var_x = 1;
% Sinal v[n]
var_v = 1e-2;
% Coeficientes do sistema
omega_0 = 0.01*pi;
omega_1 = 0.2*pi;
h0 = zeros(N,M);
h1 = zeros(N,M);
for i = 1:N
    h0_i = sin(omega_0*i);
    h1_i = sin(omega_1*i);
    for j = 1:M
        h0(i,j)=h0_i;
        h1(i,j)=h1_i;
    end
end

%Passo
mu = 0.01:0.01:0.3;
Nmu = length(mu);
% Erro médio quadrático
EMSE0_exp = zeros(Nmu,1);
EMSE1_exp = zeros(Nmu,1);
% Desvio médio quadrático
MSD0_exp = zeros(Nmu,1);
MSD1_exp = zeros(Nmu,1);
for k = 1:Nmu
    % Execuções do experimento 
    L = 100;
    % Erro médio quadrático em excesso
    EMSE0 = zeros(N,1);
    EMSE1 = zeros(N,1);
    % Desvio médio quadrático
    MSD0 = zeros(N,1);
    MSD1 = zeros(N,1);
    % Experimento
    for i = 1:L
        % Sinal x[n]
        x = sqrt(var_x)*randn(N,1);
        % Sinal v[n]
        v = sqrt(var_v)*randn(N,1);
        % LMS
        [ e0, delta_w0, ~] = LmsUnderTest(x, v, M, N, h0, mu(k), h0);
        [ e1, delta_w1, ~] = LmsUnderTest(x, v, M, N, h1, mu(k), h1);
        % Parâmetros para avaliar o aprendizado do filtro
        EMSE0 = EMSE0 + (e0.^2-var_v)/L;
        MSD0 = MSD0 + sum(delta_w0.^2,2)/L;
        EMSE1 = EMSE1 + (e1.^2-var_v)/L;
        MSD1 = MSD1 + sum(delta_w1.^2,2)/L;
    end
 
    EMSE0_exp(k) = mean(EMSE0(1500:end));
    MSD0_exp(k) = mean(MSD0(1500:end));
    EMSE1_exp(k) = mean(EMSE1(1500:end));
    MSD1_exp(k) = mean(MSD1(1500:end));
end

%%
% Plota os gráficos
figure();
plot(mu,EMSE0_exp,mu,EMSE1_exp);
title('EMSE (experimental) em função do passo \mu');
xlabel('\mu (passo do LMS)');
legend('\omega_0=0.01\pi', '\omega_0=0.2\pi');

figure();
plot(mu,MSD0_exp,mu,MSD1_exp);
title('MSD (experimental) em função do passo \mu');
xlabel('\mu (passo do LMS)');
legend('\omega_0=0.01\pi', '\omega_0=0.2\pi');