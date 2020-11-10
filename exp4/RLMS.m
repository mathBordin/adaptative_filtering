function [ e, delta_w, y ] = RLMS(x, v, h, lambda, Wo, P1)
%RLMS Algoritmo RLMS para filtros adaptativos FIR
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
Pprev = eye(M);
Pnext = eye(M);
e = zeros(N,1);
delta_w = zeros(N,M);
phi = zeros(M,1);

%%
% Saída do sistema
y = filter(h,1,x);

%%
% Sinal observado
d = y + v;

if exist('P1','var')
   Pprev = P1;
end

%%
% Iteração do algoritmo
for n = 1:N
    phi = [x(n); phi(1:M-1)];
    e(n) = d(n)-phi'*W(n,:)';
    delta_w(n,:) = Wo' - W(n,:);
    if(n~=N)
        Pnext = lambda^(-1)*Pprev - lambda^(-1)*(Pprev*phi*phi'*Pprev)/(lambda+phi'*Pprev*phi);
        W(n+1,:) = W(n,:)+(e(n)*Pnext*phi)';
        Pprev = Pnext;
    end
end

end