% CLT_Concrete_Analysis
%
% Author: MSc. Eng. Ramon Vilela
%
% Version 1.00
% 
% Updated January 28, 2024
%
% All right reserved

clc; clearvars; close all; tic

%% INPUT DATA
b   = 1.00;   % width, m
L   = 7.00;   % lenght, m
g   = 9.81;     % gravity, m/s^2
bt  = 0.14;   % lamella width, m  

% Thickness (m)
hc  = 0.10;     % concrete layer
t   = 0.00;     % spacing between materials
h1  = 0.04;     % CLT layer 1
h12 = 0.04;     % CLt layer 12
h2  = 0.04;     % CLT layer 2
h23 = 0.04;     % CLT layer 23
h3  = 0.04;     % CLT layer 3

% Concrete properties
fck = 20e6;       % characteristic strength (Pa)
alphaE = 1.2;     % factor according to NBR 6118:2023 (-)  

% Longitudinal elastic modulus of wood (Pa)
E(1) = 11000e6;    % layer 1
E(2) = 11000e6;    % layer 2
E(3) = 11000e6;    % layer 3

% Transversal elastic modulus (Pa)
G12     = 68e6;     % layer 12
G23     = 68e6;     % layer 23

% Characteristic Strength of wood (Pa)
fbk     = 24e6;     % in bending
ft0k    = 14e6;     % tension parallel
ft90k   = 0.4e6;    % tension perpendicular
fc0k    = 21e6;     % compression parallel
fc90k   = 2.5e6;    % compression perpendicular
fvk     = 4e6;      % shear
fvrk    = 1.6e6;    % rolling shear

% Mean density (kg/m^3)
rhoc    = 2400;     % concrete
rhow    = 550;      % wood

% Coefficient of creep and reduction factor
phic    = 3.4;     % concrete
phiw    = 0.6;     % wood
psi2    = 0.4;     % reduction factor according to NBR 8681:2004

% Connectors
smin    = 1.00;     % minimum spacing (m)
smax    = 1.00;     % maximum spacing (m)
nx      = 3;        % number of connectors in the x-direction of mid-plate
ny      = 5;        % number of connectors in y-direction in width b

% Loads (N/m^2)
gk      = 1e3;      % permanent
qk      = 3e3;      % live

% Variables to calculate the stiffness of the connection
m = nx;                          % number of connectors along the x-direction
n = ny;                          % number of connectors along the y-direction
% slip moduli, N/m
K = [110 80 110 80 110 
     110 80 110 80 110
     110 80 110 80 110]*10^6;
% strength of connector, N
FRs = [224 215 224 215 224
       224 215 224 215 224
       224 215 224 215 224]*0.7*10^3; 
[kser,s,x,y,Ki] = kstiffness(m, n, K, smin, smax, b, L);

% Reinforcement
d = 0.005;  % diameter of bar, m
sr = 0.10;  % spacing of reinforcement bar

%% CALCULATIONS
fid = fopen('Results.txt', 'w');    % create a results file
% Sizing
hclt  = h1+h12+h2+h23+h3;           % CLT (m)
h     = hclt + hc + t;              % CLT-Concrete (m)
nclt  = length(E);                  % quantity of layers
nlam  = floor(b/bt);                % number of lamellae

% Loads (N/m)
g1kclt = b*hclt*rhow*g;             % dead load of CLT
g1kc   = b*hc*rhoc*g;               % dead load of concrete
g2k   = gk*b;                       % permanent
q1k   = qk*b;                       % live
qSd   = 1.25*g1kclt + 1.30*g1kc + 1.35*g2k + 1.5*q1k;
qSk   = g1kclt+g1kc+g2k+q1k;

% Design bending moment (N.m)
MSd = qSd*L^2/8;  
MSk = qSk*L^2/8;

% Design shear (N)
VSd = qSd*(L/2-2*h);  

[EIclt,EAclt,gammaclt,zCG_CLT,Anet] = stiffclt(E, G12, G23, h1, h12, h2,...
    h23,h3, b, L);

% Connectors
sef   = 0.33*smin + 0.67*smax;      % spacing (m)
stot = sum(s);
V = qSd*(L/2-x);                    % shear in x-position (N)

MakeTitle(fid,'GENERAL DATA',70)
fprintf(fid,'- DIMENSIONS\n');
fprintf(fid,'Span (L): \t\t\t\t\t\t%4.0f mm\n',L*1000);
fprintf(fid,'Width (b): \t\t\t\t\t\t%4.0f mm\n',b*1000);
fprintf(fid,'Concrete thickness (hc): \t\t%4.0f mm\n',hc*1000);
fprintf(fid,'Spacing Concrete/CLT (t): \t\t%4.0f mm\n',t*1000);
fprintf(fid,'CLT thickness (hCLT): \t\t\t%4.0f mm\n',hclt*1000);
fprintf(fid,'Total thickness (h): \t\t\t%4.0f mm\n',h*1000);

fprintf(fid,'\n- LOADS\n');
fprintf(fid,'Dead load CLT (g1kCLT):\t\t\t%2.2f kN/m\n',g1kclt/1000);
fprintf(fid,'Dead load Conc. (g1kc):\t\t\t%2.2f kN/m\n',g1kc/1000);
fprintf(fid,'Permanent (g2k):\t\t\t\t%2.2f kN/m\n',g2k/1000);
fprintf(fid,'Live (qk):\t\t\t\t\t\t%2.2f kN/m\n',qk/1000);

MakeTitle(fid,'CLT STIFFNESS',70)
fprintf(fid,'Coef. of Composition (gamma1):\t%1.3f\n',gammaclt(1));
fprintf(fid,'Coef. of Composition (gamma2):\t%1.3f\n',gammaclt(2));
fprintf(fid,'Coef. of Composition (gamma3):\t%1.3f\n',gammaclt(3));
fprintf(fid,'Centroide (zCG):\t\t\t\t%3.0f mm\n',zCG_CLT*1000);
fprintf(fid,'Bending Stiffness (EIclt):\t\t%1.3fe6 Nm^2\n',EIclt*10^-6);

MakeTitle(fid,'SPACING OF CONNECTORS',70)
fprintf(fid,'Minimum spacing (smin): \t\t%2.0f mm\n',smin*1000);
fprintf(fid,'Maximum spacing (smin): \t\t%2.0f mm\n',smax*1000);
fprintf(fid,'smax <= 4smin: \t\t\t\t\t');
if smax <= 4*smin
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'smax <= 1000 mm: \t\t\t\t');
if smax <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'Effective spacing (sef): \t\t%2.0f mm\n',sef*1000);

%% ULTIMATE LIMIT STATE (ULS)
% Connector
FRds1 = FRs./1.4;           % Design strength of connector, N
k = 2/3*kser;               % distributed stiffness of connector, N/m^2

% Concrete
if fck <= 50e6                      % mean tension strength (Pa)
    fctm = 0.3*(fck*10^-6)^(2/3)*10^6;
else
    fctm = 2.12*log(1 + 0.11*fck*10^-6)*10^6;
end
fctkinf = 0.7*fctm;         % inferior characteristic strength (Pa)
fctd = fctkinf/1.4;         % design tension strength (Pa)
fcd = fck/1.4;              % design compression strength (Pa)
Eci = alphaE*5600*sqrt(fck*10^-6)*10^6; 
                            % initial Young's modulus (Pa)
alphai = 0.8 + 0.2*fck/80e6;  
                            % calculation factor
Ecs = alphai*Eci;           % secant Young's modulus (Pa)
nbar = b/sr;                % number of bars
As1 = nbar*pi*d^2/(4);      % area of reinforcement, m^2
tauRd = 0.25*fctd;          % design shear strength (Pa)
sigmacp = 0;                % normal stress (Pa)

MakeTitle(fid,'MECHANICAL PROPERTIES',70)
fprintf(fid,'- Compression strength of concrete:\n');
fprintf(fid,'Characteristic (fck): \t\t\t%2.2f MPa\n',fck/10^6);
fprintf(fid,'Design (fcd): \t\t\t\t\t%2.2f MPa\n',fcd/10^6);
fprintf(fid,'- Tensile strength of concrete:\n');
fprintf(fid,'Medium (fctm): \t\t\t\t\t%2.2f MPa\n',fctm/10^6);
fprintf(fid,'Inferior (fctk,inf): \t\t\t%2.2f MPa\n',fctkinf/10^6);
fprintf(fid,'Design (fctd): \t\t\t\t\t%2.2f MPa\n',fctd/10^6);
fprintf(fid,'- Shear strength of concrete:\n');
fprintf(fid,'Design (tauRd): \t\t\t\t%2.2f MPa\n',tauRd/10^6);
fprintf(fid,'\n- Elasticit modulus of concrete:\n');
fprintf(fid,'Initial (Eci): \t\t\t\t\t%2.2f GPa\n',Eci/10^9);
fprintf(fid,'Secant (Ecs): \t\t\t\t\t%2.2f GPa\n',Ecs/10^9);

% Madeira
kmod = 0.7*1*0.95;          % modification factor
ksys  = (nlam + 34)/35;     % modification system factor
ft0d = ksys*kmod*ft0k/1.4;  % design parallel tension strength (Pa)
fbd = ksys*kmod*fbk/1.4;    % design bending strength (Pa)
fvrd = kmod*fvrk/1.8;       % rolling shear strength (Pa)

% Calcs of (EI)ef
[EIef,hcef,a,gamma,EA,EI,it] = EIeffective(E,Ecs,G12,G23,fctd,k,hc,t,h1,...
    h12,h2,h23,h3,b,L,MSd);
fprintf(fid,'Iterations: \t\t\t\t\t%.0f\n',it);

% Discretization along the thickness
z(1) = 0;
z(2) = 0;
z(3) = h1;
z(4) = h1;
z(5) = h1+h12;
z(6) = h1+h12;
z(7) = h1+h12+h2;
z(8) = h1+h12+h2;
z(9) = h1+h12+h2+h23;
z(10) = h1+h12+h2+h23;
z(11) = hclt;
z(12) = hclt;
z(13) = h-hcef;
z(14) = h-hcef;
z(15) = h;
z(16) = h;

% Intermediate stress (Pa)
sigma_csup = -MSd*Ecs/EIef*(gamma(1)*a(1)+0.5*hcef);	% top concrete
sigma_cinf = -MSd*Ecs/EIef*(gamma(1)*a(1)-0.5*hcef);	% underside conc.
sigma_cltsup = MSd*E(3)/EIef*(gamma(2)*a(2)-0.5*hclt);	% top CLT
sigma_cltinf = MSd*E(1)/EIef*(gamma(2)*a(2)+0.5*hclt);	% underside CLT
sigma = zeros(1,16);
for i=2:11
    sigma(i) = (sigma_cltsup-sigma_cltinf)*z(i)/hclt + sigma_cltinf;
end
sigma(1) = 0;
sigma(4) = 0;
sigma(5) = 0;
sigma(8) = 0;
sigma(9) = 0;
sigma(12) = 0;
sigma(13) = 0;
sigma(14) = sigma_cinf;
sigma(15) = sigma_csup;
sigma(16) = 0;

zLNclt = -sigma_cltinf*hclt/(sigma_cltsup-sigma_cltinf);
zLNc = h-hcef-sigma_cinf*hcef/(sigma_csup-sigma_cinf);

% Strain (m/m)
strain = sigma./E(1);
strain(14) = sigma_cinf/Ecs;
strain(15) = sigma_csup/Ecs;

% Figure of stress along the section
figure(3)
hold on
plot(sigma*10^-6,z*10^3)
plot([min(sigma) max(sigma)]*10^-6,[zLNclt zLNclt]*10^3,'--',...
    'Color',[.6,.3,0])
if zLNc <= hclt
else
    plot([min(sigma) max(sigma)]*10^-6,[zLNc zLNc]*10^3,'--',...
        'Color',[.5,.5,.5])
end
plot([0 0],[0 h*10^3],'k')
box on
set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
ylim([0 h*10^3])
xlabel('Normal stress (MPa)')
ylabel('CLT-Concrete thickness (mm)')
legend('Stress', 'NA CLT', 'NA Concrete','EdgeColor','none')
title('NORMAL STRESS AT MIDSPAN - ULS')
saveas(gcf,'StressResults.png')

%% RESISTANCE BENDING MOMENT
nxPoints = ceil(L/.2);
nyPoints = ceil(b/.1);
xi = linspace(0, L/2, nxPoints);
yi = linspace(-b/2, b/2, nyPoints);
for i=1:length(xi)
    for j=1:length(yi)
        M(i,j) = -qSd*xi(i)^2/2 + qSd*L^2/8;
    end
end

% Bendin moment resistance of the Concrete
MRdc = fcd*EIef/(Ecs*(0.5*hc+gamma(1)*a(1)));   
etabc = M./MRdc;
etabcmax = max(etabc(:));
tit = sprintf('\\eta_{b,c} = \t%.2f', etabcmax);
PlotDesignRatio(etabc, b, L, xi, yi, tit);
note{1} = tit;

% Bending moment resistence of the CLT
MRdclt = EIef/(0.5*E(1)*hclt/fbd + a(2)*EA(2)/(ft0d*Anet));
etabclt = M./MRdclt;
etabcltmax = max(etabclt(:));
tit = sprintf('\\eta_{b,CLT} = \t%.2f', etabcltmax);
PlotDesignRatio(etabclt, b, L, xi, yi, tit);
note{2} = tit;

MRd = min([MRdc,MRdclt]);            % bending moment resistence
etab = MSd/MRd;                      % design factor
etabc = MSd/MRdc;                    % design factor - concrete
etabclt = MSd/MRdclt;                % design factor - CLT

fprintf(fid,'Concrete thickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac)\t%.3f\n',gamma(1));
fprintf(fid,'\n- BENDING MOMENT\n');
fprintf(fid,'Design bending moment (MSd): \t%.2f kN.m\n',MSd/1000);
fprintf(fid,'Resistant bending CLT (MRdCLT): %.2f kN.m\n',MRdclt/1000);
fprintf(fid,'Resistant bending Conc. (MRdc): %.2f kN.m\n',MRdc/1000);
fprintf(fid,'Design factor (etab): \t\t\t%.2f %%\n',etab*100);

%% RESISTANCE SHEAR EFFORT
Fs = zeros(m,n);
Fsi = zeros(m,1);
for i=1:m
    for j=1:n
        Fsi(i) = gamma(1)*EA(1)*a(1)*s(1)*qSd*x(i)/EIef;
        Fs(i,j) = K(i,j)/Ki(i)*Fsi(i);
        %Fs(i,j) = K(i,j)/Ki(i)*gamma(1)*EA(1)*a(1)*s(1)*qSd*x(i)/EIef;
    end
end
etavs = Fs./FRds1;
etavsmax = max(etavs(:));
tit = sprintf('\\eta_{v,s} = \t%.2f', etavsmax);
PlotDesignRatio(etavs, b, L, x, y, tit);
note{3} = tit;

% Maximum shear
VRds = VSd/etavsmax;                 % in function of the connector

% Design factor of shear in the CLT
for i=1:length(xi)
    for j=1:length(yi)
        V(i,j) = qSd*xi(i);
    end
end
VRclt = EI(2)*b*fvrd/(E(1)*b*h1*(h1/2+h12+h2/2));
VRdclt = EIef/(EI(2)+0.5*EA(2)*a(2)*hclt)*VRclt;
etavclt = V./VRdclt;
etavcltmax = max(etavclt(:));
tit = sprintf('\\eta_{v,CLT} = \t%.2f', etavcltmax);
PlotDesignRatio(etavclt, b, L, xi, yi, tit);
note{4} = tit;

% Design factor to shear in Concrete
rho1 = As1/(b*hcef);                   % reinforcement rate 
VRdc1 = (tauRd*(1.2+40*rho1)+0.15*sigmacp)*b*hcef;
VRdc = EIef/(EI(1)+0.5*gamma(1)*EA(1)*a(1)*(2*hc-hcef+t))*VRdc1;
etavc = V./VRdc;
etavcmax = max(etavc(:));
tit = sprintf('\\eta_{v,c} = \t%.2f', etavcmax);
PlotDesignRatio(etavc, b, L, xi, yi, tit);
note{5} = tit;
                                
VRd = min([VRds,VRdclt,VRdc]);       % shear resistance of concrete

% Design factor
etav = VSd/VRd;                      % general
etavs = VSd/VRds;                   % design factor - connector
etavclt = VSd/VRdclt;                % design factor - CLT
etavc = VSd/VRdc;                    % design factor - Concrete

fprintf(fid,'\n- SHEAR EFFORT\n');
fprintf(fid,'Design shear effort (VSd): \t\t%.2f kN\n',VSd/1000);
fprintf(fid,'Resistant shear connec. (VRds): %.2f kN\n',VRds/1000);
fprintf(fid,'Resistant shear CLT (VRdCLT):\t%.2f kN\n',VRdclt/1000);
fprintf(fid,'Resistant shear Conc. (VRdc):\t%.2f kN\n',VRdc/1000);
fprintf(fid,'Design factor (etav): \t\t\t%.2f %%\n',etav*100);

%% SERVICE LIMIT STATE - INSTANTANEOUS DISPLACEMENT
% Calcs of (EI)ef
[EIef,hcef,a,gamma,EA,EI,it] = ...
    EIeffective(E,Ecs,G12,G23,fctd,kser,...
    hc,t,h1,h12,h2,h23,h3,b,L,MSk);

% Loads
p = [g1kclt    
    g1kc
    g2k
    q1k];

% Displacements
winst = sum(5.*p.*L^4./(384*EIef)); % instantaneous displacement (m)
winstlim = L/500;                  	% limit of instantaneous  
                                  	% displacement (m)
etawinst = winst/winstlim;        	% design factor

% Print numerical results
fprintf(fid,'\n:::::::::::::SERVICE LIMIT STATE::::::::::::::\n');
fprintf(fid,'Number of iterations: \t\t\t%.0f\n',it);
fprintf(fid,'Concrete thickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac):\t%.3f\n',gamma(1));
fprintf(fid,'\n- INSTANTANEOUS DISPLACEMENT\n');
fprintf(fid,'Displacement (winst): \t\t\t%.1f mm\n',winst*1000);
fprintf(fid,'Limit (winstlim): \t\t\t\t%.1f mm\n',winstlim*1000);
fprintf(fid,'Design factor (etawinst): \t\t%.2f %%\n',etawinst*100);

% Design ratio to instantaneous displacement
for i=1:length(xi)
    for j=1:length(yi)
        ii = length(xi)+1-i;
        w(ii,j) = sum(p)*(10*L^4-43*L^3*xi(i)+48*L^2*xi(i)^2-8*xi(i)^4)...
            /(768*EIef)-5*sum(p)*L^4/(384*EIef); 
    end
end
etawinst = abs(w)./winstlim;
etawinstmax = max(etawinst(:));
tit = sprintf('\\eta_{w,inst} = \t%.2f', etawinstmax);
PlotDesignRatio(etawinst, b, L, xi, yi, tit);
note{6} = tit;

%% Service Limit State - Vibration
ml = rhow*b*hclt+rhoc*b*hc+g2k/g;   % linear mass (kg/m)
f1 = pi/(2*L^2)*sqrt(EIef/ml);      % natural frequency (Hz)
w1kN = 1000*L^3/(48*EIef);          % displacement duo to 1 kN load
w1kNlim = f1^1.43/39000;            % limit displacement
etaf = w1kN/w1kNlim;                % design factor

fprintf(fid,'\n- EXCESSIVE VIBRATION\n');
fprintf(fid,'Natural frequency (f1): \t\t%.2f Hz\n',f1);
fprintf(fid,'Displacement 1kN (w1kN): \t\t%.2f mm\n',w1kN*1000);
fprintf(fid,'Disp. limit (w1kNlim): \t\t\t%.2f mm\n',w1kNlim*1000);
fprintf(fid,'Design factor (etaf): \t\t\t%.2f %%\n',etaf*100);

% Design ratio to natural frequency
for i=1:length(xi)
    for j=1:length(yi)
        %ii = length(xi)+1-i;
        w(i,j) = 1000*(L^3-6*L*xi(i)^2+4*xi(i)^3)/(48*EIef); 
        %z = length(xi)+1-i;
    end
end
etaw1kN = abs(w)./w1kNlim;
etaw1kNmax = max(etaw1kN(:));
tit = sprintf('\\eta_{vib} = \t%.2f', etaw1kNmax);
PlotDesignRatio(etaw1kN, b, L, xi, yi, tit);
note{7} = tit;

%% Service Limit State - Final displacement
% Calcs of (EI)ef
k = kser./(2.*(1+phiw));            % creep connector stiffness
Ec = Ecs/(1+phic);                  % creep Young's modulus of concrete
E = E./(1+phiw);                    % creep elastic modulus in bending 
                                    % of wood layers
G12 = G12./(1+phiw);                % creep transversal elastic modulus                                   
                                    % of 12 wood layer   
G23 = G23./(1+phiw);                % creep transversal elastic modulus                                   
                                    % of 23 wood layers
[EIef,hcef,a,gamma,EA,EI,it] = ...
    EIeffective(E,Ec,G12,G23,fctm,k,...
    hc,t,h1,h12,h2,h23,h3,b,L,MSk);

p(4) = p(4)*psi2;                   % reduced live load (N/m)
wfin = sum(5.*p.*L^4./(384*EIef));  % final displacement (m)
wfinlim = L/300;                    % limit displacement (m)
etawfin = wfin/wfinlim;             % design factor

fprintf(fid,'\n- FINAL DISPLACEMENT\n');
fprintf(fid,'Number of iterations: \t\t\t%.0f\n',it);
fprintf(fid,'Concrete thickness (hcef): \t\t%.2f mm\n',hcef*1000);
fprintf(fid,'Eff. bending stiffness (EI)ef:');
fprintf(fid,'\t%2.3fe6 N.m^2\n',EIef*10^-6);
fprintf(fid,'Coef. of Composition (gammac):\t%.3f\n',gamma(1));
fprintf(fid,'Displacement (wfin): \t\t\t%.1f mm\n',wfin*1000);
fprintf(fid,'Limit (wfinlim): \t\t\t\t%.1f mm\n',wfinlim*1000);
fprintf(fid,'Design factor (etawfin): \t\t%.2f %%\n',etawfin*100);

% Design ratio to instantaneous displacement
for i=1:length(xi)
    for j=1:length(yi)
        ii = length(xi)+1-i;
        w(ii,j) = sum(p)*(10*L^4-43*L^3*xi(i)+48*L^2*xi(i)^2-8*xi(i)^4)...
            /(768*EIef)-5*sum(p)*L^4/(384*EIef); 
    end
end
etawfin = abs(w)./wfinlim;
etawfinmax = max(etawfin(:));
tit = sprintf('\\eta_{w,fin} = \t%.2f', etawfinmax);
PlotDesignRatio(etawfin, b, L, xi, yi, tit);
note{8} = tit;

%% PLOT TABLE OF RESULTS
% ULS-Bending moment
fprintf(fid,'\nRESUME');
fprintf(fid,'\n---------------------------------');
fprintf(fid,'-------------------------------\n');
fprintf(fid,'CRITERION \t\t\t\t | DF \t\t | GOVERN \t | ACCEPTANCE\n');
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');
% Concrete
fprintf(fid,'ULS - Bending moment\t | %.2f %%\t | ',etabc*100);
fprintf(fid,'Concrete  | ');
if etabc <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% CLT
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etabclt*100);
fprintf(fid,'CLT \t\t | ');
if etabclt <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% ULS-Shear Effort
fprintf(fid,'ULS - Shear effort\t\t | %.2f %%\t | ',etavs*100);
% Connector
fprintf(fid,'Connector | ');
if etavs <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% Concrete
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etavc*100);
fprintf(fid,'Concrete  | ');
if etavc <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
% CLT
fprintf(fid,'\t\t\t\t\t\t | %.2f %%\t | ',etavclt*100);
fprintf(fid,'CLT \t\t | ');
if etavclt <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% SLS-Instantaneous displacement 
fprintf(fid,'SLS - Inst. Displacement | %.2f %%\t | ',max(etawinst(:))*100);
fprintf(fid,'--- \t\t | ');
if max(etawinst(:)) <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% SLS-Final displacement
fprintf(fid,'SLS - Final displacement | %.2f %%\t | ',max(etawfin(:))*100);
fprintf(fid,'--- \t\t | ');
if max(etawfin(:)) <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

% SLS-Vibration
fprintf(fid,'SLS - Vibration \t\t | %.2f %%\t | ',etaf*100);
fprintf(fid,'--- \t\t | ');
if etaf <= 1
    fprintf(fid,'OK!\n');
else
    fprintf(fid,'NOT ACCEPTED!\n');
end
fprintf(fid,'---------------------------------');
fprintf(fid,'-------------------------------\n');

fclose(fid);                    % close fid

figure(1)
hold on
etamax = [etabcmax, etabcltmax, etavcmax, etavcltmax, etavsmax,...
    etaw1kNmax, etawfinmax, etawinstmax];
[maxEta, maxIndex] = max(etamax);

% Nomes correspondentes aos etas
etaNames = {'b,c', 'b,CLT', 'v,c', 'v,CLT', 'v,s', 'fn', 'w,fin',...
    'w,inst'};
maxEtaName = etaNames{maxIndex};

% Configuração do título
title(sprintf('MAXIMUM DESIGN RATIO \\eta_{%s} = %.2f', maxEtaName, maxEta));
txt = text(L/4, 0, 1.5*maxEta, note, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
    'FontSize', 14, 'Color', 'Black','BackgroundColor', 'white');
set(txt, 'EdgeColor', 'Black', 'BackgroundColor', [1, 1, 1, 0.5]);

clearvars                       % Clear variables

%% OPEN RESULTS
open Results.txt
fprintf('\nTotal time:\t%.2f s\n',toc)
