function [k,s,x,y,Ki] = kstiffness(m, n, K, smin, smax, b, L)
% Calculates stiffness of connectors
% Input data
%   m     = number of connector positions along x
%   n     = number of connector rows along y
%   K     = matrix of slip moduli of connectors at position ij, N/m
%   smin  = minimum spacing between connectors in x-axis, m
%   smax  = maximum spacing between connectors in x-axis, m
%   b     = width of the specimen, m
%   L     = total length of the specimen, m
%
% Output results
%   k     = distributed stiffness of conectors, N/m^2
%   s     = vector of spacing, m
%   x     = vector of position of connectors in x-axis, m
%   y     = vector of position of connectors in y-axis, m
%   Ki    = vector of equivalente slip modulus to each position in x-axis,
%           N/m

% Initializes vector s
s = zeros(1,m);
x = zeros(1,m);
y = zeros(1,n);
for i = 1:m
    s(i) = (smin-smax)/(1-m)*(i-m)+smax;
    if i==1
        x(i) = (L-s(i))/2;
    else
        x(i) = x(i-1) - (s(i-1)+s(i))/2;
    end
    for j=1:n
        y(j) = b/n*(j-(1+n)/2);
    end
end

% Initializes variable k
k = 0;
Ki = zeros(m,1);
% Calculates double sum based on the provided equation
for i = 1:m
    for j = 1:n
        k = k + K(i, j) / (s(i) * m);
    end
    Ki(i) = sum(K(i,:));
end
% % Surface of slip moduli
% [X, Y] = meshgrid(x, y);
% figure;
% surf(X, Y, K', 'EdgeColor', 'none');
% xlabel('Mid length (m)');
% ylabel('Width (m)');
% zlabel('Slip modulus');
% title('Slip moduli in function of surface');
% xlim([0 L/2]);
% ylim([-b/2 b/2]);
% zlim([0 max(K(:))]);
% grid on;
% box on;
end