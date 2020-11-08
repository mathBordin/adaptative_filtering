function [ e, delta_w, W ] = LmsUnderTest(x, v, M, N, h, mu, Wo)
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
delta_w = zeros(N,M);
x_aux = zeros(M,1);
y_est = zeros(N,1);
% Saída do sistema
y = variableFirFilter(h,x);
% Sinal observado
d = y + v;

szWo = size(Wo);
Wo_t=zeros(N,M);
if(szWo(1) == 1)
    for i = 1:N
        Wo_t(i,:) = Wo;
    end
else
    Wo_t = Wo;
end

        
% Iteração do algoritmo
for n = 1:N
    x_aux = [x(n); x_aux(1:M-1)];
    y_est(n) = W(n,:)*x_aux;
    e(n) = d(n)-y_est(n);
    delta_w(n,:) = Wo_t(n,:) - W(n,:);
    if(n~=N)
        W(n+1,:) = W(n,:)+mu*e(n)*x_aux';
    end
end

end
