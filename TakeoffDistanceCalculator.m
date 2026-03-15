%Variables

%Zero-lift drag coefficient
Cd0 = 0.025;
%Wing area (m^2)
S = 16.2;
%Aspect ratio
Ar = 7.32;
%Number of simultaneous equations
n = 10;
%Values of theta
theta = linspace(10*pi()/180, 80*pi()/180, n);
%Wing span(m)
b = 11;
%Lift curve slope (1/rad)
a0 = 2*pi();
%Root chord (m)
c_root = 1.63;
%Tip chord (m)
c_tip = 1.13;
%Span location where tapering begins (m)
taper_loc = 2.68;
%Aircraft mass (kg)
mass = 1043;
%Acceleration due to gravity (m/s^2)
g = 9.81;
%Air density (kg/m^3)
rho = 1.225;
%Time step (s)
dt = 0.05;
%Maximum lift coefficient of airfoil
CLmax = 1.3;
%Stall velocity (m/s)
V_s = sqrt((2*mass*g)/(rho*CLmax*S));
%Propeller efficiency
eta = 0.8;
%Engine power (W)
P = 119312;
%Propeller radius (m)
Rp = 0.95;
%Propeller disk area (m^2)
Ap = pi*Rp^2;
%Empty velocity and distance vectors
V = [0,0];
S_t = [0,0];
%Matrices to input calculated values into
B = zeros(n, 1);
X = zeros(n, n);
nu = zeros(n, 1);
c = zeros(n, 1);

%Lifting line theory
for i = 1:n
    L = 0.5 * b * cos(theta);
    if L(i) > taper_loc
        c(i) = (c_root - c_tip)*(0.5*b*(1-cos(theta(i))))/(0.5*b - taper_loc) + c_tip;
    else
        c(i) = c_root;
    end
    nu(i) = (c(i)*a0)/(4*b) ;
    alpha_0 = -pi()/180;
    B(i) = nu(i) * sin(theta(1, i)) * -(alpha_0);
    for j = 1:n
        X(i, j) = sin((2*j - 1)*theta(1, i))*(sin(theta(1, i)) + (2*j - 1)*nu(i));
    end
end

solutions = rhs\lhs;

C_l = solutions(1)*pi()*Ar;

delta = 0;
for i = 2:n
    delta = delta + (2*i - 1)*(solutions(i)/solutions(1))^2;
end

e = 1/(1 + delta);
C_d = Cd0 + (C_l^2/(pi()*e*Ar));

i = 1
K1 = (2*eta*P)/(rho*Ap);
K4 = -(rho*S/(2*mass))*(Cd0 + (C_l^2)/(pi*e*Ar));

while V(i) < 1.2*V_s
    K2 = real(-V(i)/3);
    Ve = real((K2^3 + K1 + sqrt(2 * K2^3 * K1 + K1^2))^(1/3) + (K2^3 + K1 - sqrt(2 * K2^3 * K1 + K1^2))^(1/3) + K2);
    T = real(0.5 * rho * Ap * (Ve + 3*K2) * (Ve - 3*K2));
    K3 = real((T/mass));
    i = i + 1;
    V(i) = 1/(K4*dt) + sqrt(1/(K4^2 * dt^2) - V(i - 1)^2 - (2*K3)/K4 - (2*V(i - 1))/(K4*dt));
    S_t(i) = S_t(i - 1) + (V(i) + V(i - 1))*dt/2;
end

plot(S_t, V)
S_takeoff = max(real(S_t))

%Empty acceleration, velocity and distance vectors
a = [0,0];
V = [0,0];
S_t = [0,0];

i = 1;
while V(i) < 1.2*V_s
    K2 = real(-V(i)/3);
    Ve = real((K2^3 + K1 + sqrt(2 * K2^3 * K1 + K1^2))^(1/3) + (K2^3 + K1 - sqrt(2 * K2^3 * K1 + K1^2))^(1/3) + K2);
    T = real(0.5 * rho * Ap * (Ve + 3*K2) * (Ve - 3*K2));
    K3 = real((T/mass));
    i = i + 1;
    a(i) = K3 + K4*V(i - 1)^2;
    V(i) = V(i - 1) + a(i - 1)*dt;
    S_t(i) = S_t(i - 1) + V(i - 1)*dt;
end

plot(S_t, V)