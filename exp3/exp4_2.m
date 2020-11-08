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
a = 0.9;
Hb_in=(1-a^2);
Ha_in=[1 -a];
var_x = 1;
% Sinal v[n]
var_v = 1e-2;
% Coeficientes do sistema
var_q = 1e-4;
H = [1 0.5];
% Coeficientes do filtro adaptativo
W = zeros(N,M);
% Matrix de autocorrelação
R_phi = [(1-a^2) a*(1-a^2); a*(1-a^2) (1-a^2)]; 
% Passo que minimiza o EMSE
mu_o = sqrt((M*var_q)/(var_v*trace(R_phi)));

%%
% Calculo teórico do EMSE e do MSD em regime
mu = 0.01:0.01:0.3;
Nmu = length(mu);
EMSE_teo = (mu*var_v*trace(R_phi)/2)+(M*var_q./(2*mu));
EMSE_o = (mu_o*var_v*trace(R_phi)/2)+(M*var_q./(2*mu_o));
MSD_teo = (mu*M*var_v/2)+(var_q*trace(inv(R_phi))./(2*mu));
MSD_o = (mu_o*M*var_v/2)+(var_q*trace(inv(R_phi))./(2*mu_o));

%% 
% Estimação experimental do EMSE e do MSD em regime

% Erro médio quadrático
EMSE_exp = zeros(Nmu,1);
% Desvio médio quadrático
MSD_exp = zeros(Nmu,1);
for k = 1:Nmu
    % Execuções do experimento 
    L = 100;
    % Erro médio quadrático em excesso
    EMSE = zeros(N,1);
    % Desvio médio quadrático
    MSD = zeros(N,1);
    % Experimento
    for i = 1:L
        % Sinal x[n]
        x = sqrt(var_x)*randn(N,1);
        x = filter(Hb_in,Ha_in,x);
        % Sinal v[n]
        v = sqrt(var_v)*randn(N,1);
        % Sistema H
        h = createFilterPertubation( H, N, var_q );
        % LMS
        [ e, delta_w, ~] = LmsUnderTest(x, v, M, N, h, mu(k), h);
        % Parâmetros para avaliar o aprendizado do filtro
        EMSE = EMSE + (e.^2-var_v)/L;
        MSD = MSD + sum(delta_w.^2,2)/L;
    end
    EMSE_exp(k) = mean(EMSE(500:end));
    MSD_exp(k) = mean(MSD(500:end));
end

%%
% Plota os gráficos
figure();
plot(mu,EMSE_teo,mu,EMSE_exp);
hold on;
plot(mu_o,EMSE_o,'rx');
legend('Teórico','Experimental', 'Ótimo');
title('EMSE em função do passo \mu');

figure();
plot(mu,MSD_teo,mu,MSD_exp);
hold on;
plot(mu_o,MSD_o,'rx');
legend('Teórico','Experimental', 'Ótimo');
title('MSD em função do passo \mu');