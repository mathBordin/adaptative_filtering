function [ W, e, y_est, emse, msd ] = LMS_under_test(x, v, M, N, H, mu, Wot)
%LMS Algoritmo LMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x é o sinal de referência
% v é o sinal de interesse (em caso real, é desconhecido)
% M é a ordem do filtro adaptativo
% N é o tamanho dos sinais
% Coeficientes do sistema FIR de interesse
% mu é o passo do algoritmo
% W é o vetor com a evolução temporal do algoritmo
% e é o vetor de erro
% y_est é a saída do filtro adaptativo

% Declaração de variáveis iniciais
W = zeros(N,M);  
e = zeros(N,1);
emse = zeros(N,1);
msd = zeros(N,1);
x_aux = zeros(M,1);
y_est = zeros(N,1);
% Saída do sistema
y = filter(H, 1, x);
% Sinal observado
d = y + v;

% Iteração do algoritmo
for n = 1:N
    x_aux = [x(n); x_aux(1:M-1)];
    y_est(n) = W(n,:)*x_aux;
    e(n) = d(n)-y_est(n);
    delta_w = Wot - W(n,:);
    emse(n) = (delta_w*x_aux)^2;
    msd(n) = delta_w*delta_w';
    if(n~=N)
        W(n+1,:) = W(n,:)+mu*e(n)*x_aux';
    end
end

end
