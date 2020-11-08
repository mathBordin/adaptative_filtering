function tan = mean_tan(b, e, arr)
%mean_tan Fun��o para calcular a tangente m�dia entre os pontos b e e 
%   Matheus Bordin Gomes
% b � o ponto inicial
% e � o ponto final
% arr � o vetor de dados
% tan � a tangente m�dia

N = floor((e-b)/2-1);
tan_arr = zeros(N,1);

for i = 0:N-1
    tan_arr(i+1) = (arr(b+i*2+2) - arr(b+i*2))/2;
end

tan = sum(tan_arr)/N;

end
