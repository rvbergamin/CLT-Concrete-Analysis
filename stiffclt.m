function [EIclt,EAclt,gammaclt,zCG_CLT,Anet] = stiffclt(E, G12, G23, h1, h12, h2, h23, h3, b, L)

% This function calculates the stifness properties of CLT based on
% Input data
%   E     = vector of elastic bending modulus of lamellae, N/m^2
%   G12   = shear modulus of layer 12, N/m^2
%   G23   = shear modulus of layer 23, N/m^2
%   hi    = tickness of layer i = 1, 2 or 3, m
%   hij   = tickness of layer ij = 12 or 23, m
%   b     = width, m
%   L     = span, m
%
% Output results
%   EIclt = bending stiffness of CLT, N.m^2
%   EAclt = axial stiffness of CLT, N
%   gammaclt = coeficiente of composition to longitudinal layers
%   zCG_CLT = central of gravity of CLT, m

nclt  = length(E);                  % quantity of layers

% CLT areas (m^2)
A(1)      = b*h1;                   % Area da camada 1
A(2)      = b*h2;                   % Area da camada 2
A(3)      = b*h3;                   % Area da camada 3

% CLT moment of inertia (m^4)
I(1)      = b*h1^3/12;              % Area da camada 1
I(2)      = b*h2^3/12;              % Area da camada 2
I(3)      = b*h3^3/12;              % Area da camada 3

% Composition coefficient (-)
gammaclt(1)   = 1/(1+pi^2*E(1)*A(1)*h12/(L^2*b*G12));
gammaclt(2)   = 1;
gammaclt(3)   = 1/(1+pi^2*E(3)*A(3)*h23/(L^2*b*G23));

% Position CG's layers in relation of base(m)
zclt(1) = h1/2;
zclt(2) = h1+h12+h2/2;
zclt(3) = h1+h12+h2+h23+h3/2;

yEAz = 0;
yEA = 0;
for i=1:nclt
    yEAz = gammaclt(i)*E(i)/E(1)*A(i)*zclt(i) + yEAz;
    yEA = gammaclt(i)*E(i)/E(1)*A(i) + yEA;
end
zCG_CLT = yEAz/yEA;

aclt = zeros(1,nclt);
for i=1:nclt
    aclt(i) = zclt(i)-zCG_CLT;
end

% (EI)clt
EIclt = 0;
for i=1:nclt
    EIclt = E(i)*(I(i)+gammaclt(i)*A(1)*aclt(i)^2)+ EIclt;
end

% (EA)clt
EAclt = 0;
for i=1:nclt
    EAclt = E(i)*A(i)+EAclt;
end
Anet = sum(A);             % Net area (m^2)
end