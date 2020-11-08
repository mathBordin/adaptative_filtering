function y = variableFirFilter(h,x)
% Lê tamanho do filtro e do sinal
[Ns,~] = size(x);
[Nf,M] = size(h);
% Utiliza o menor tamanho de amostras para realizar a filtragem
if Ns<Nf
    N = Ns;
else 
    N = Nf;
end
y = zeros(N,1);
x_aux = zeros(M,1);
% Filtragem do sinal
for n=1:N    
    x_aux = [x(n); x_aux(1:M-1)];
    y(n) = h(n,:)*x_aux;
end

end

