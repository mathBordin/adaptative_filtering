function [ w, e, delta_w_0, delta_w_1 ] = LMS(N, x, v, H_b, mi)

% Vetor que conterá a evolução dos coeficientes do filtro adaptativo 
w = zeros(N,2);  
% Saída do sistema
y = filter(H_b, 1, x);
% Sinal observado
d = y + v;
% Estimativa de y
y_est = zeros(N,1);
% Erro
e = zeros(N,1);
% Primeira iteração do algoritmo adaptativo com condição inicial nula
y_est(1) =  w(1, :)*[x(1) 0]';
e(1) = d(1) - y_est(1);
w(2,:) = w(1,:) + mi*e(1)*[x(1) 0];
% Iterações seguintes do algoritmo adaptativo
for i=2:N
    y_est(i) =  w(i, :)*[x(i) x(i-1)]';
    e(i) = d(i) - y_est(i);
    if(i ~= N)
        w(i+1,:) = w(i,:) + mi*e(i)*[x(i) x(i-1)];
    end
end
% Variação de cada um dos coeficientes do filtro adaptativo em cada passo
delta_w_0 = H_b(1)-w(:,1);
delta_w_1 = H_b(2)-w(:,2);

end
