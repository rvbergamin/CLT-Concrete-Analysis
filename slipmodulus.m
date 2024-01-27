% Valores das variáveis de exemplo
m = 3;  % número de conectores ao longo de x
n = 4;  % número de conectores ao longo de y
K = [80.31 110.24 80.31 110.24
     80.31 110.24 80.31 110.24
     80.31 110.24 80.31 110.24];  % módulo de deslizamento, N/m
s = [(0.3+0.5)/2
     (0.5+0.5)/2
     (0.5+0.7)/2];

% Inicializando a variável k
k = 0;

% Calculando a soma dupla com base na equação fornecida
for i = 1:m
    for j = 1:n
        k = k + K(i, j) / (s(i) * m);
    end
end

% Exibindo o resultado
disp(['O valor de k é: ' num2str(k)]);
