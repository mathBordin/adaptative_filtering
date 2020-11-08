%%
% Experiência 1: Curvas de aprendizado
% PSI3531 - Prof. Vítor
% Matheus Bordin Gomes - 9838028
clear;clc;close all;
%#ok<*NOPTS>

%%
% Exercício 1 - (omega_0 = 0.5*pi)

N = 2000; % Número de amostras dos sinais
n = 1:N;  % Vetor de amostras 
omega_0 = 0.5*pi;
omega_1 = 0.01*pi; 
A = 0.001;
% Sinal interferente
Theta = rand()*2*pi;
x = cos(omega_0*n+Theta); 
% Sinal de interesse
Psi = rand()*2*pi;
v = A*cos(omega_1*n+Psi); 
% Sistema
H_b = [1; -0.5];  
% Passo do algoritmo
mu = 0.01; 

% Execução do algoritmo
[ w, e, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
delta_w_0 = H_b(1) - w(:,1);
delta_w_1 = H_b(2) - w(:,2);

% Plota os gráficos
figure(1);
plot(n, delta_w_0, n, delta_w_1);
title('Variação dos coeficientes do filtro adaptativo (\omega_0 = 0.5*pi)');
xlabel('Passos do algoritmo');
legend('\Delta w_0', '\Delta w_1');

figure(2);
plot(n, e);
title('Erro da estimação com o filtro adaptativo (\omega_0 = 0.5*pi)');
xlabel('Passos do algoritmo');
ylabel('Erro');

figure(3);
plot(n, w(:,1), n, w(:,2), n, H_b(1)*ones(size(n)), n, H_b(2)*ones(size(n)));
title('Evolução dos coeficientes do filtro adaptativo (omega_0 = 0.5*pi)');
xlabel('Passos do algoritmo');
legend('w_0','w_1', 'H_0', 'H_1');

% Matriz de autocorrelação e seus autovalores
disp('Resultados para omega_0 = 0.5*pi');
R_phi = [0.5 0.5*cos(omega_0); 0.5*cos(omega_0) 0.5]
eigenvalues = eig(R_phi);
fprintf('Autovalores de R_phi: %.5f %.5f \n', eigenvalues(1), eigenvalues(2));

% Taxa de convergência
taxa_w0_exp = mean_tan(1,1400,log(abs(delta_w_0)));
taxa_w1_exp = mean_tan(1,1400,log(abs(delta_w_1)));
taxa_e2_exp = mean_tan(1,800,log(e.^2)); 
taxa_teorica = eig(eye(2)-mu*R_phi);
fprintf('Taxa de convergência experimental de delta_w0: %.5f \n', taxa_w0_exp);
fprintf('Taxa de convergência experimental de delta_w1: %.5f \n', taxa_w1_exp);
fprintf('Taxa de convergência experimental de e^2: %.5f \n', taxa_e2_exp);
fprintf('Valor de eig(I-mu*R_phi): %.5f %.5f\n\n\n', taxa_teorica(1), taxa_teorica(2) );

%%
% Exercício 1 - (omega_0 = 0.1*pi)

N = 30000; % Número de amostras dos sinais
n = 1:N;  % Vetor de amostras 

omega_0 = 0.1*pi;
% Sinal interferente
Theta = rand()*2*pi;
x = cos(omega_0*n+Theta); 
% Sinal de interesse
Psi = rand()*2*pi;
v = A*cos(omega_1*n+Psi); 

% Execução do algoritmo
[ w, e, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
delta_w_0 = H_b(1) - w(:,1);
delta_w_1 = H_b(2) - w(:,2);

% Plota os gráficos
figure(4);
plot(n, delta_w_0, n, delta_w_1);
title('Variação dos coeficientes do filtro adaptativo (omega_0 = 0.1*pi)');
xlabel('Passos do algoritmo');
legend('\Delta w_0', '\Delta w_1');

figure(5);
plot(n, e);
title('Erro da estimação com o filtro adaptativo (omega_0 = 0.1*pi)');
xlabel('Passos do algoritmo');
ylabel('Erro');

figure(6);
plot(n, w(:,1), n, w(:,2), n, H_b(1)*ones(size(n)), n, H_b(2)*ones(size(n)));
title('Evolução dos coeficientes do filtro adaptativo (omega_0 = 0.1*pi)');
xlabel('Passos do algoritmo');
legend('w_0','w_1', 'H_0', 'H_1');

% Matriz de autocorrelação e seus autovalores
disp('Resultados para omega_0 = 0.1*pi');
R_phi = [0.5 0.5*cos(omega_0); 0.5*cos(omega_0) 0.5]
eigenvalues = eig(R_phi);
fprintf('Autovalores de R_phi: %.5f %.5f \n', eigenvalues(1), eigenvalues(2));

% Taxa de convergência
taxa_w0_exp = mean_tan(5000,20000,log(abs(delta_w_0)));
taxa_w1_exp = mean_tan(5000,20000,log(abs(delta_w_1)));
taxa_e2_exp = mean_tan(5000,10000,log(e.^2)); 
taxa_teorica = eig(eye(2)-mu*R_phi);

fprintf('Taxa de convergência experimental de delta_w0: %.5f \n', taxa_w0_exp);
fprintf('Taxa de convergência experimental de delta_w1: %.5f \n', taxa_w1_exp);
fprintf('Taxa de convergência experimental de e^2: %.5f \n', taxa_e2_exp);
fprintf('Valor de eig(I-mu*R_phi): %.5f %.5f\n\n\n', taxa_teorica(1), taxa_teorica(2) );
figure(13)
plot(n,log(abs(delta_w_0)))

%%
% Exercício 2 - Curva de aprendizado (omega_0 = 0.1*pi)')

% Curva de aprendizado experimental
P = 30; % Número de realizações
E_delta_w_0_exp = zeros(N,1); % Determinação experimental de dalta_w0
E_delta_w_1_exp = zeros(N,1); % Determinação experimental de dalta_w1
R_phi = [0.5 0.5*cos(omega_0); 0.5*cos(omega_0) 0.5];
for i = 1:P
    % Sinal interferente
    Theta = rand()*2*pi;
    x = cos(omega_0*n+Theta); 
    % Sinal de interesse
    Psi = rand()*2*pi;
    v = A*cos(omega_1*n+Psi); 
    %Execução do algoritmo
    [ w, ~, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
    delta_w_0 = H_b(1) - w(:,1);
    delta_w_1 = H_b(2) - w(:,2);
    E_delta_w_0_exp = E_delta_w_0_exp+delta_w_0;
    E_delta_w_1_exp = E_delta_w_1_exp+delta_w_1;
end
E_delta_w_0_exp = E_delta_w_0_exp/P; 
E_delta_w_1_exp = E_delta_w_1_exp/P;

% Curva de aprendizado teórica
E_delta_w_teo = zeros(N,2);
E_delta_w_teo(1,:) = [E_delta_w_0_exp(1) E_delta_w_1_exp(1)];
for i = 2:N
    E_delta_w_teo(i,:) = (eye(2)-mu*R_phi)*E_delta_w_teo(i-1,:)';
end

% Plota as curvas
figure(7);
plot(n, E_delta_w_0_exp, n, E_delta_w_teo(:,1));
title('Curvas de aprendizado de w0 teórica x experimental (omega_0 = 0.1*pi)');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w0');

figure(8);
plot(n, E_delta_w_1_exp, n, E_delta_w_teo(:,2));
title('Curvas de aprendizado de w1 teórica x experimental (omega_0 = 0.1*pi)');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w1');

%%
% Exercício 2 - Curva de aprendizado (omega_0 = 0.5*pi)')

N = 2000; % Número de amostras dos sinais
n = 1:N;  % Vetor de amostras 

% Curva de aprendizado experimental
P = 30; % Número de realizações
E_delta_w_0_exp = zeros(N,1); % Determinação experimental de dalta_w0
E_delta_w_1_exp = zeros(N,1); % Determinação experimental de dalta_w1
omega_0 = 0.5*pi;
R_phi = [0.5 0.5*cos(omega_0); 0.5*cos(omega_0) 0.5];
for i = 1:P
    % Sinal interferente
    Theta = rand()*2*pi;
    x = cos(omega_0*n+Theta); 
    % Sinal de interesse
    Psi = rand()*2*pi;
    v = A*cos(omega_1*n+Psi); 
    %Execução do algoritmo
    [ w, ~, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
    delta_w_0 = H_b(1) - w(:,1);
    delta_w_1 = H_b(2) - w(:,2);
    E_delta_w_0_exp = E_delta_w_0_exp+delta_w_0;
    E_delta_w_1_exp = E_delta_w_1_exp+delta_w_1;
end
E_delta_w_0_exp = E_delta_w_0_exp/P; 
E_delta_w_1_exp = E_delta_w_1_exp/P;

% Curva de aprendizado teórica
E_delta_w_teo = zeros(N,2);
E_delta_w_teo(1,:) = [E_delta_w_0_exp(1) E_delta_w_1_exp(1)];
for i = 2:N
    E_delta_w_teo(i,:) = (eye(2)-mu*R_phi)*E_delta_w_teo(i-1,:)';
end

% Plota as curvas
figure(9);
plot(n, E_delta_w_0_exp, n, E_delta_w_teo(:,1));
title('Curvas de aprendizado de w0 teórica x experimental (omega_0 = 0.5*pi)');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w0');

figure(10);
plot(n, E_delta_w_1_exp, n, E_delta_w_teo(:,2));
title('Curvas de aprendizado de w1 teórica x experimental (omega_0 = 0.5*pi)');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w1');

%%
% Exercício 3 - Caso em que os sinais são ruídos brancos

% Sinal interferente
x = randn(N,1); 
% Sinal de interesse
var_n = 0.01;
v = sqrt(var_n)*randn(N,1); 

% Execução do algoritmo
[ w, e, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
delta_w_0 = H_b(1) - w(:,1);
delta_w_1 = H_b(2) - w(:,2);

% Plota os gráficos
figure(11);
plot(n, delta_w_0, n, delta_w_1);
title('Variação dos coeficientes do filtro adaptativo - sinais WGN');
xlabel('Passos do algoritmo');
legend('\Delta w_0', '\Delta w_1');

figure(12);
plot(n, e);
title('Erro da estimação com o filtro adaptativo - sinais WGN');
xlabel('Passos do algoritmo');
ylabel('Erro');

figure(13);
plot(n, w(:,1), n, w(:,2), n, H_b(1)*ones(size(n)), n, H_b(2)*ones(size(n)));
title('Evolução dos coeficientes do filtro adaptativo - sinais WGN');
xlabel('Passos do algoritmo');
legend('w_0','w_1', 'H_0', 'H_1');

% Matriz de autocorrelação e seus autovalores
disp('Resultados para caso com sinais WGN');
r = xcorr(x,1,'biased');
R_phi = [r(2) r(1); r(1) r(2)]
eigenvalues = eig(R_phi);
fprintf('Autovalores de R_phi: %.5f %.5f \n', eigenvalues(1), eigenvalues(2));

% Taxa de convergência
taxa_w0_exp = mean_tan(1,300,log(abs(delta_w_0)));
taxa_w1_exp = mean_tan(1,300,log(abs(delta_w_1)));
taxa_e2_exp = mean_tan(1,300,log(e.^2)); 
taxa_teorica = eig(eye(2)-mu*R_phi);
fprintf('Taxa de convergência experimental de delta_w0: %.5f \n', taxa_w0_exp);
fprintf('Taxa de convergência experimental de delta_w1: %.5f \n', taxa_w1_exp);
fprintf('Taxa de convergência experimental de e^2: %.5f \n', taxa_e2_exp);
fprintf('Valor de eig(I-mu*R_phi): %.5f %.5f\n\n', taxa_teorica(1), taxa_teorica(2) );

% Curva de aprendizado experimental
P = 30; % Número de realizações
E_delta_w_0_exp = zeros(N,1); % Determinação experimental de dalta_w0
E_delta_w_1_exp = zeros(N,1); % Determinação experimental de dalta_w1
for i = 1:P
    % Sinal interferente
    x = randn(N,1); 
    % Sinal de interesse
    var_n = 0.01;
    v = sqrt(var_n)*randn(N,1);  
    %Execução do algoritmo
    [ w, ~, ~] = LMS_under_test(x, v, 2, N, H_b, mu); 
    delta_w_0 = H_b(1) - w(:,1);
    delta_w_1 = H_b(2) - w(:,2);
    E_delta_w_0_exp = E_delta_w_0_exp+delta_w_0;
    E_delta_w_1_exp = E_delta_w_1_exp+delta_w_1;
end
E_delta_w_0_exp = E_delta_w_0_exp/P; 
E_delta_w_1_exp = E_delta_w_1_exp/P;

% Curva de aprendizado teórica
E_delta_w_teo = zeros(N,2);
E_delta_w_teo(1,:) = [E_delta_w_0_exp(1) E_delta_w_1_exp(1)];
for i = 2:N
    E_delta_w_teo(i,:) = (eye(2)-mu*R_phi)*E_delta_w_teo(i-1,:)';
end

% Plota as curvas
figure(14);
plot(n, E_delta_w_0_exp, n, E_delta_w_teo(:,1));
title('Curvas de aprendizado de w0 teórica x experimental');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w0');

figure(15);
plot(n, E_delta_w_1_exp, n, E_delta_w_teo(:,2));
title('Curvas de aprendizado de w1 teórica x experimental');
legend('Experimental', 'Teórica');
xlabel('Passos do algoritmo');
xlabel('Delta w1');

%%
% Cálculo das variâncias
var_delta_w0_exp = var(delta_w_0(1200:end));
var_delta_w1_exp = var(delta_w_1(1200:end));
var_delta_w_teo = 0.5*mu*var_n;
var_delta_e0_exp = var(e(1200:end)-v(1200:end));
fprintf('var(delta_w0) experimental: %.5f \n', var_delta_w0_exp);
fprintf('var(delta_w1) experimental: %.5f \n', var_delta_w1_exp);
fprintf('var(delta_w) teórico: %.5f \n', var_delta_w_teo);
fprintf('var(e-v) experimental: %.5f \n', var_delta_e0_exp);