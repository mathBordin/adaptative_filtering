function [ e, delta_w, W ] = LmsUnderTest(x, v, M, N, h, mu, Wo)
%LMS Algoritmo LMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x � o sinal de refer�ncia
% v � o sinal de interesse (em caso real, � desconhecido)
% M � a ordem do filtro adaptativo
% N � o tamanho dos sinais
% Coeficientes do sistema FIR de interesse
% mu � o passo do algoritmo
% W � o vetor com a evolu��o temporal do algoritmo
% e � o vetor de erro
% y_est � a sa�da do filtro adaptativo

% Declara��o de vari�veis iniciais
W = zeros(N,M);  
e = zeros(N,1);
delta_w = zeros(N,M);
x_aux = zeros(M,1);
y_est = zeros(N,1);
% Sa�da do sistema
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

        
% Itera��o do algoritmo
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
