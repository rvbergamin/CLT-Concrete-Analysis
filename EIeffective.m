function [EIef,hcef,a,gamma,EA,EI,it,tef,sigma] = EIeffective(E,Ec,G12,...
    G23,fct,k,hc,t,h1,h12,h2,h23,h3,b,L,M)

% This function calculates the stifness properties of CLT based on
% Input data
%   E     = vector of elastic bending modulus of lamellae, N/m^2
%   Ec    = elastic modulus of concrete, N/m^2
%   G12   = shear modulus of layer 12, N/m^2
%   G23   = shear modulus of layer 23, N/m^2
%   fct   = tensile strength of concrete, N/m^2
%   
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

hclt  = h1+h12+h2+h23+h3;           % CLT (m)
h = hc + t + hclt;                  % total tickness, m 
hcef=hc;                            % effective thickness (m)
sigma_cinf = fct;                   % stress inderside of concrete (Pa)
it=0;                               % iteration counter
EA(2) = 0;                          % axial stiffness of CLT (N)
nclt  = length(E);                  % quantity of layers
sigma = 0;

[EIclt,EAclt] = stiffclt(E, G12, G23, h1, h12, h2, h23, h3, b, L);

EI(2) = EIclt;
EA(2) = EAclt;

% Loop to reduce effective concrete thickness
while sigma_cinf >= fct
    if sigma_cinf <= fct || it==0
        hcef=hc;
    else
        hcef=hcef*(1-.005);
    end
    
    EA(1) = Ec*b*hcef;         % Axial stiffness of concrete, N
    
    % Composition coefficient (-)
    %gamma(1) = 1/(1+pi^2*EA(1)*sef/(K*L^2)); % Concreto
    gamma(1) = 1/(1+pi^2*EA(1)/(k*L^2)); % Concreto
    gamma(2) = 1;                            % CLT
    
    % Distance between CG (m)
    at = h-(hcef+hclt)/2;
    a(1) = gamma(2)*EA(2)/(gamma(1)*EA(1)+gamma(2)*EA(2))*at;
    a(2) = gamma(1)*EA(1)/(gamma(1)*EA(1)+gamma(2)*EA(2))*at;
    
    % Effective bending stiffness (N.m^2)
    EI(1) = Ec*b*hcef^3/12;
    EIef = 0;
    for i=1:2
        EIef = EI(i)+ gamma(i)*EA(i)*a(i)^2 + EIef;
    end
    
    % Stress underside conc.(Pa)
    sigma_cinf = -M*Ec/EIef*(gamma(1)*a(1)-0.5*hcef);
    sigma(it+1) = sigma_cinf;
    it=it+1;
end
tef = t+ hc-hcef;
end