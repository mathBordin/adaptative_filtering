function [ e, delta_w, y ] = NLMS(x, v, h, mu, eps, Wo)
%NLMS Algoritmo NLMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x é o sinal de referência
% v é o sinal de interesse (em caso real, é desconhecido)
% h são os coeficientes do sistema FIR de interesse (podem ou não variar no tempo)
% mu é o passo do algoritmo
% eps é o fator que evita a normalização do passo por valores muito pequenos
% W0 são os coeficientes ótimos do filtro (podem ou não variar no tempo)
% e é o vetor de erro
% delta_w é a diferença entre o Wo e o W
% W são os coeficientes do filtro adaptativo

%%
% Define o tamanho dos vetores utilizados

% Tamanho dos sinais
if length(x) > length(v)
    N = length(x);
else
    N = length(v);
end
% Comprimento do filtro
M = length(h);

%%
% Declaração de variáveis iniciais
W = zeros(N,M);  
e = zeros(N,1);
delta_w = zeros(N,M);
x_aux = zeros(M,1);
y_est = zeros(N,1);

%%
% Saída do sistema
y = filter(h,1,x);

%%
% Sinal observado
d = y + v;
        
% Iteração do algoritmo
for n = 1:N
    x_aux = [x(n); x_aux(1:M-1)];
    y_est(n) = W(n,:)*x_aux;
    e(n) = d(n)-y_est(n);
    delta_w(n,:) = Wo' - W(n,:);
    if(n~=N)
        W(n+1,:) = W(n,:)+(mu/(x_aux'*x_aux+eps))*e(n)*x_aux';
    end
end

end
