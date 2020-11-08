function tan = mean_tan(b, e, arr)
%mean_tan Função para calcular a tangente média entre os pontos b e e 
%   Matheus Bordin Gomes
% b é o ponto inicial
% e é o ponto final
% arr é o vetor de dados
% tan é a tangente média

N = floor((e-b)/2-1);
tan_arr = zeros(N,1);

for i = 0:N-1
    tan_arr(i+1) = (arr(b+i*2+2) - arr(b+i*2))/2;
end

tan = sum(tan_arr)/N;

end
