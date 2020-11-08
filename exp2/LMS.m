function [W,e,y]=LMS(x, d, M, N, mu)
%LMS Algoritmo LMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x é o sinal de referência
% d é o sinal observado
% M é a ordem do filtro adaptativo
% N é o tamanho dos sinais
% mu é o passo do algoritmo
% W é o vetor com a evolução temporal do algoritmo
% e é o vetor de erro
% y é a saída do filtro adaptativo

% Declaração de variáveis iniciais
x_aux=zeros(M,1);
W=zeros(N,M); 
y=zeros(N,1);
e=zeros(N,1);

% Iteração do algoritmo
for n = 1:N
    x_aux = [x(n); x_aux(1:M-1)];
    y(n) = W(n,:)*x_aux;
    e(n) = d(n)-y(n);
    if(n~=N)
        W(n+1,:) = W(n,:)+mu*e(n)*x_aux';
    end
end

end