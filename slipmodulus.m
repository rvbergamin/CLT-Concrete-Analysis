% Valores das vari�veis de exemplo
m = 3;  % n�mero de conectores ao longo de x
n = 4;  % n�mero de conectores ao longo de y
K = [80.31 110.24 80.31 110.24
     80.31 110.24 80.31 110.24
     80.31 110.24 80.31 110.24];  % m�dulo de deslizamento, N/m
s = [(0.3+0.5)/2
     (0.5+0.5)/2
     (0.5+0.7)/2];

% Inicializando a vari�vel k
k = 0;

% Calculando a soma dupla com base na equa��o fornecida
for i = 1:m
    for j = 1:n
        k = k + K(i, j) / (s(i) * m);
    end
end

% Exibindo o resultado
disp(['O valor de k �: ' num2str(k)]);
