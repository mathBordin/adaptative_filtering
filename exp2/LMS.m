function [W,e,y]=LMS(x, d, M, N, mu)
%LMS Algoritmo LMS para filtros adaptativos FIR
%   Matheus Bordin Gomes
% x � o sinal de refer�ncia
% d � o sinal observado
% M � a ordem do filtro adaptativo
% N � o tamanho dos sinais
% mu � o passo do algoritmo
% W � o vetor com a evolu��o temporal do algoritmo
% e � o vetor de erro
% y � a sa�da do filtro adaptativo

% Declara��o de vari�veis iniciais
x_aux=zeros(M,1);
W=zeros(N,M); 
y=zeros(N,1);
e=zeros(N,1);

% Itera��o do algoritmo
for n = 1:N
    x_aux = [x(n); x_aux(1:M-1)];
    y(n) = W(n,:)*x_aux;
    e(n) = d(n)-y(n);
    if(n~=N)
        W(n+1,:) = W(n,:)+mu*e(n)*x_aux';
    end
end

end