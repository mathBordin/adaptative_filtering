function [ e, delta_w, y, run_time ] = RLS(x, v, h, lambda, Wo, P1)
%RLS Algoritmo RLS para filtros adaptativos FIR
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
P = eye(M);
e = zeros(N,1);
delta_w = zeros(N,M);
phi = zeros(M,1);

%%
% Sa�da do sistema
y = filter(h,1,x);

%%
% Sinal observado
d = y + v;

if exist('P1','var')
   P = P1;
end

%%
% Itera��o do algoritmo
tic
for n = 1:N
    phi = [x(n); phi(1:M-1)];
    e(n) = d(n)-phi'*W(n,:)';
    delta_w(n,:) = Wo' - W(n,:);
    if(n~=N)
        P_phi = P*phi;
        P_phi_n = P_phi/sqrt(lambda*(lambda+phi'*P_phi));
        P = lambda^(-1)*P/lambda - (P_phi_n*P_phi_n');
        W(n+1,:) = W(n,:)+(e(n)*P*phi)';
    end
end
run_time = toc;

end