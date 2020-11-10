function [ e, delta_w, y ] = NLMS(x, v, h, mu, eps, Wo)
%NLMS Algoritmo NLMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x � o sinal de refer�ncia
% v � o sinal de interesse (em caso real, � desconhecido)
% h s�o os coeficientes do sistema FIR de interesse (podem ou n�o variar no tempo)
% mu � o passo do algoritmo
% eps � o fator que evita a normaliza��o do passo por valores muito pequenos
% W0 s�o os coeficientes �timos do filtro (podem ou n�o variar no tempo)
% e � o vetor de erro
% delta_w � a diferen�a entre o Wo e o W
% W s�o os coeficientes do filtro adaptativo

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
% Declara��o de vari�veis iniciais
W = zeros(N,M);  
e = zeros(N,1);
delta_w = zeros(N,M);
x_aux = zeros(M,1);
y_est = zeros(N,1);

%%
% Sa�da do sistema
y = filter(h,1,x);

%%
% Sinal observado
d = y + v;
        
% Itera��o do algoritmo
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
