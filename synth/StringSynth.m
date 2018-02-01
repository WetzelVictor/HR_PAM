close all;
clearvars;
clc;

% String parameters
f1 = 196; % G3
f1 = f1*(2^(-5/12));
d_s = 4.4e-4; % d'Addario Steel G (4th) string
% d_s = 8.5e-4; % Alliance Savarez Nylon G (4th) string
dens_s = mean([7.75e3,8.05e3]); % Density of steel in kg/m^3
% dens_s = 1.15e3; % Density of nylon in kg/m^3
rho_s = dens_s*(0.25*pi*d_s^2); % Linear mass density of string
L_s = 0.655; % Vibrating length of string
E_s = 1.6153e11; % Steel Young modulus from measurement analysis
% E_s = 1.1916e9; % Nylon Young modulus from measurement analysis
I_s = (pi*(d_s^4))/64; % Cylinder moment of inertia for string
T_s = (rho_s*(2*L_s*f1)^2)-(((pi^2)*E_s*I_s)/(L_s^2)); % Required tension
kap = sqrt(T_s/(E_s*I_s)); % String stiffness constant
Z = sqrt(T_s*rho_s); % String driving-point impedance
B = ((pi^2)*E_s*I_s)/(T_s*(L_s^2));

% Pick parameters
L_p = 0.03; % Pick length
Weqv_p = 0.0141; % Width of equivalent rectangle (area equal to pick)
h_p = 9.6e-4; % Pick thickness (bar height)
I_p = (Weqv_p*(h_p^3))/12; % Second moment of area for pick as bar
a = 1.5e-3; % Half-width of pick point
E_p = 4e9; % High end of range for Young modulus of Nylon
% alpha = (2*Z*(L_p^3))/(3*E_p*I_p);


% Coupling parameters
mu = 0.1; % Static friction coefficient, Nylon to steel
% mu = mean([0.15,0.25]); % Static friction coefficient, Nylon to Nylon
exAng = struct('down',75,'up',15);
theta = deg2rad(exAng.up); % Excitation angle
for l = 1:10
    alpha = (2*Z*(L_p^3))/(3*E_p*I_p);
    P = ((4*Z*L_p)/(3*alpha))*mu; % Excitation force
    d = (P+0.1795)/(1.396e3);
    E_p = (P*(L_p^3))/(3*d*I_p);
end
% P = ((4*Z*L_p)/(3*alpha))*mu; % Excitation force
P_v = P*cos(theta); % Vertical polarisation (w.r.t. table)
P_h = P*sin(theta); % Horizontal polarisation
exPos = struct('P1',0.56,'P2',0.518,'P3',0.468);
x0 = exPos.P2; % Excitation position
x0_m = x0-a;
x0_p = x0+a;

% Global parameters and constants
Nx = 1000; % Number of discrete spatial points to work with
delx = L_s/(Nx-1);
x = 0:delx:L_s; % Spatial coordinate vector
rho_a = 1.184; % Air density at ~300K
eta = 1.85e-5; % Air viscosity coefficient at 300K
Qth = 3e4; % Thermo-elastic dissipation factor
del_ve = 6e-3; % Visco-elastic dissipation factor, steel
% del_ve = 4e-3; % Visco-elastic dissipation factor, Nylon

% Initial string deformation
M_L = ((x0_p/L_s)-(sinh(kap*x0_p)/sinh(kap*L_s)))*((kap/tanh(kap*L_s))-(1/L_s));
M_L = M_L-((((L_s-x0_m)/L_s)-(sinh(kap*(L_s-x0_m))/sinh(kap*L_s)))*((1/L_s)-(kap/sinh(kap*L_s))));
M_L = M_L/((((kap/tanh(kap*L_s))-(1/L_s))^2)-((1/L_s)-((kap/sinh(kap*L_s)))^2));
M0 = (((L_s-x0_m)/L_s)-(sinh(kap*(L_s-x0_m))/sinh(kap*L_s)))*((kap/tanh(kap*L_s))-(1/L_s));
M0 = M0-(((x0_p/L_s)-(sinh(kap*x0_p)/sinh(kap*L_s)))*((1/L_s)-(kap/sinh(kap*L_s))));
M0 = M0/((((kap/tanh(kap*L_s))-(1/L_s))^2)-((1/L_s)-((kap/sinh(kap*L_s)))^2));
y0_v = zeros(Nx,1);
y0_h = y0_v;
Nx0_m = floor(x0_m/delx)+1;
Nx0_p = floor(x0_p/delx)+1;
vmTerm1 = ((P_v*(L_s-x0_m))/(T_s*L_s))*x(1:Nx0_m);
vmTerm2 = ((P_v*sinh(kap*(L_s-x0_m)))/(kap*T_s*sinh(kap*L_s)))*sinh(kap*x(1:Nx0_m));
% vmTerm3 = ((P_v*M_L)/(T_s*L_s))*x(1:Nx0_m);
% vmTerm4 = ((P_v*M_L)/(T_s*sinh(kap*L_s)))*sinh(kap*x(1:Nx0_m));
% vmTerm5 = ((P_v*M0)/(T_s*L_s))*(L_s-x(1:Nx0_m));
% vmTerm6 = ((P_v*M0)/(T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(1:Nx0_m)));
y0_v(1:Nx0_m) = vmTerm1-vmTerm2;
vpTerm1 = ((P_v*x0_p)/(T_s*L_s))*(L_s-x(Nx0_p:Nx));
vpTerm2 = ((P_v*sinh(kap*x0_p))/(kap*T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(Nx0_p:Nx)));
% vpTerm3 = ((P_v*M_L)/(T_s*L_s))*x(Nx0_p:Nx);
% vpTerm4 = ((P_v*M_L)/(T_s*sinh(kap*L_s)))*sinh(kap*x(Nx0_p:Nx));
% vpTerm5 = ((P_v*M0)/(T_s*L_s))*(L_s-x(Nx0_p:Nx));
% vpTerm6 = ((P_v*M0)/(T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(Nx0_p:Nx)));
y0_v(Nx0_p:Nx) = vpTerm1-vpTerm2;
y0_v(Nx0_m:Nx0_p)=((y0_v(Nx0_p)-y0_v(Nx0_m))/(x(Nx0_p)-x(Nx0_m)))*(x(Nx0_m:Nx0_p)-x(Nx0_m))+y0_v(Nx0_m);
y0_v(isnan(y0_v)) = 0;
hmTerm1 = ((P_h*(L_s-x0_m))/(T_s*L_s))*x(1:Nx0_m);
hmTerm2 = ((P_h*sinh(kap*(L_s-x0_m)))/(kap*T_s*sinh(kap*L_s)))*sinh(kap*x(1:Nx0_m));
hmTerm3 = ((P_h*M_L)/(T_s*L_s))*x(1:Nx0_m);
hmTerm4 = ((P_h*M_L)/(T_s*sinh(kap*L_s)))*sinh(kap*x(1:Nx0_m));
hmTerm5 = ((P_h*M0)/(T_s*L_s))*(L_s-x(1:Nx0_m));
hmTerm6 = ((P_h*M0)/(T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(1:Nx0_m)));
y0_h(1:Nx0_m) = hmTerm1-hmTerm2-hmTerm3+hmTerm4-hmTerm5+hmTerm6;
hpTerm1 = ((P_h*x0_p)/(T_s*L_s))*(L_s-x(Nx0_p:Nx));
hpTerm2 = ((P_h*sinh(kap*x0_p))/(kap*T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(Nx0_p:Nx)));
hpTerm3 = ((P_h*M_L)/(T_s*L_s))*x(Nx0_p:Nx);
hpTerm4 = ((P_h*M_L)/(T_s*sinh(kap*L_s)))*sinh(kap*x(Nx0_p:Nx));
hpTerm5 = ((P_h*M0)/(T_s*L_s))*(L_s-x(Nx0_p:Nx));
hpTerm6 = ((P_h*M0)/(T_s*sinh(kap*L_s)))*sinh(kap*(L_s-x(Nx0_p:Nx)));
y0_h(Nx0_p:Nx) = hpTerm1-hpTerm2-hpTerm3+hpTerm4-hpTerm5+hpTerm6;
y0_h(Nx0_m:Nx0_p)=((y0_h(Nx0_p)-y0_h(Nx0_m))/(x(Nx0_p)-x(Nx0_m)))*(x(Nx0_m:Nx0_p)-x(Nx0_m))+y0_h(Nx0_m);
y0_h(isnan(y0_h)) = 0;

% Modal synthesis
Nm = 25; % Number of modes to include in synthesis
T = 7; % Time period over which to synthesize
fs = 44100; % Sampling frequency
Ts = 1/fs; % Sampling period
t = 0:Ts:T; % Time vector
Nt = length(t);
y_v = zeros(Nx,Nt);
y_h = zeros(Nx,Nt);
fn = zeros(Nm,1);
for n = 1:Nm
    k = (n*pi)/L_s;
    Y_v = sin(k*x);
    Y_h = sin(k*x)+((k/kap)*(cosh(kap*x)-cos(k*x)-sinh(kap*x)));
    Y_h(isnan(Y_h)) = 0;
    b_v = (1/trapz(abs(Y_v).^2))*trapz(y0_v*Y_v);
    b_h = (1/trapz(abs(Y_h).^2))*trapz(y0_h*Y_h);
    f = (n/(2*L_s))*sqrt((T_s/rho_s)*(1+(B*(n^2))));
    fn(n) = f;
    omega = 2*pi*f;
    M = (d_s/4)*sqrt(omega/eta);
    R = (pi^2)*rho_a*f*((d_s^2)/2)*((sqrt(2)/M)+(1/(2*(M^2))));
    Qinv = (R/(omega*rho_s))+(((omega^2)*rho_s*E_s*I_s*del_ve)/(T_s^2))+(1/Qth);
    y_t = exp(-pi*f*Qinv*t).*cos(omega*t);
    y_v = y_v+b_v*(Y_v.')*y_t;
    y_h = y_h+b_h*(Y_h.')*y_t;
end
y = y_v+y_h;

yx0_m = y(Nx0_m,:);
yx0_m = yx0_m/max(abs(yx0_m));

% snr = 43.8746; % Signal-to-noise ratio in dB for Nylon string measurements
snr = 46.8754; % Signal-to-noise ratio in dB for steel string measurements
G = 10^(-snr/20);
g = G*((2*rand(1,Nt))-1);

if (rms(yx0_m((Nt-round(fs/f1)):Nt))>=G)
    while (rms(yx0_m((Nt-round(fs/f1)):Nt))>=G)
        T = T+1;
        t = 0:Ts:T; % Time vector
        Nt = length(t);
        y_v = zeros(Nx,Nt);
        y_h = zeros(Nx,Nt);
        fn = zeros(Nm,1);
        for n = 1:Nm
            k = (n*pi)/L_s;
            Y_v = sin(k*x);
            Y_h = sin(k*x)+((k/kap)*(cosh(kap*x)-cos(k*x)-sinh(kap*x)));
            Y_h(isnan(Y_h)) = 0;
            b_v = (1/trapz(abs(Y_v).^2))*trapz(y0_v*Y_v);
            b_h = (1/trapz(abs(Y_h).^2))*trapz(y0_h*Y_h);
            f = (n/(2*L_s))*sqrt((T_s/rho_s)*(1+(B*(n^2))));
            fn(n) = f;
            omega = 2*pi*f;
            M = (d_s/4)*sqrt(omega/eta);
            R = (pi^2)*rho_a*f*((d_s^2)/2)*((sqrt(2)/M)+(1/(2*(M^2))));
            Qinv = (R/(omega*rho_s))+(((omega^2)*rho_s*E_s*I_s*del_ve)/(T_s^2))+(1/Qth);
            y_t = exp(-pi*f*Qinv*t).*cos(omega*t);
            y_v = y_v+b_v*(Y_v.')*y_t;
            y_h = y_h+b_h*(Y_h.')*y_t;
        end
        y = y_v+y_h;
        yx0_m = y(Nx0_m,:);
        yx0_m = yx0_m/max(abs(yx0_m));
        g = G*((2*rand(1,Nt))-1);
    end
end

% save('Table1/Steel/PickG_UpP1','fn','t','y_v','y_h','fs');
clearvars -except yx0_m g t fs;

sig = yx0_m+g;
sig = sig/max(abs(sig));

% figure;
% plot(t,sig,'LineWidth',2);
% xlabel('Temps (s)');
% ylabel('Amplitude');
% ax = gca;
% ax.LineWidth = 2;
% ax.FontSize = 16;
% fig = gcf;
% pos = fig.PaperPosition;
% fig.PaperSize = [pos(3),pos(4)];
% print('Table1/Steel/PickG_UpP1','-dpdf');

audiowrite('Instru/Steel/D.wav',sig,fs);