% 1. Original Transfer Function
%{
num=6.68;
den=[114.89 91.44 32];
t=0:0.01:20;
step(num,den,t)
%}

% 2. PID Controller
% 2.A. Initialization of Inputs
rise = 0.25; % Rise time
overshoot = 0.75; % Overshoot percentage
settling = 0.35; % Settling time

% Extra treatments for input variables to align
% with stepinfo results from 2.B. (modifiers chosen via trial-and-error)
rise = rise - 0.13;
overshoot = overshoot - 11;
settling = settling + 0.691;

% 2.B. Initialization of PID variables
% Initialization of x and n values for the equation estimation
% via multivariable linear regression (with 3 input variables & 1 output)
nval = 12;
x1val = transpose([5.0514 2.5095 1.2303 0.1002 43.1582 26.1001 40.6735 52.8634 0.3272 0.2204 0.1328 0.1721]);
x2val = transpose([0 0 0 4.0692 0.9322 0.7747 0.3513 0.0918 9.2012 1.8565 1.9999 0]);
x3val = transpose([18.8623 12.4115 9.9226 0.6693 55.0804 33.8965 53.6264 71.9802 5 3.0792 1.6125 0.2968]);

% Equation estimation for Kd
syms x;
yval1 = transpose([1 15.5 25 325 30 15 20 25 101 205 320 215]);
eqny = triinputeqn(nval,yval1,x1val,x2val,x3val);
Kd1 = eqny(1) + (eqny(2)*rise) + (eqny(3)*overshoot) + (eqny(4)*settling); % Multivariable linear regression
Kd2 = double(subs((-7836.7347*(x^3))+(12149.2711*(x^2))-(6109.2128*x)+1064.6099, x, settling)); % Eqn2 (Excel)
Kd3 = double(subs((4.5605*(x^2))-(85.544*x)+373.4686, x, overshoot)); % Eqn3 (Excel)
% Average of Kd1-3 to fit more in the PID characteristic relationship
Kd = max((Kd1 + Kd2 + Kd3) / 3, 0);

% Equation estimation for Kp
yval1 = transpose([5 14.7 27 650 0.67 0.83 0.23 0.18 100 46 53 197]);
eqny = triinputeqn(nval,yval1,x1val,x2val,x3val);
Kp1 = eqny(1) + (eqny(2)*rise) + (eqny(3)*overshoot) + (eqny(4)*settling); % Multivariable linear regression
Kp2 = double(subs((-440000*(x^3))+(330285.7143*(x^2))-(84114.2857*x)+7593, x, rise)); % Eqn1 (Excel)
Kp3 = double(subs((-0.1581*(x^2))+(20.9542*x)-196.8111, x, overshoot)); % Eqn4 (Excel)
% Average of Kp1-3 to fit more in the PID characteristic relationship
Kp = max((Kp1 + Kp2 + Kp3) / 3, 0);

% Equation estimation for Ki
yval1 = transpose([1.3 3 5 105 0.27 0.43 0.26 0.2 128 35 80.55 3 ]);
eqny = triinputeqn(nval,yval1,x1val,x2val,x3val);
Ki1 = eqny(1) + (eqny(2)*rise) + (eqny(3)*overshoot) + (eqny(4)*settling); % Multivariable linear regression
Ki2 = double(subs((-440000*(x^3))+(330285.7143*(x^2))-(84114.2857*x)+7593, x, rise)); % Eqn1 (Excel)
Ki3 = double(subs((-7836.7347*(x^3))+(12149.2711*(x^2))-(6109.2128*x)+1064.6099, x, settling)); % Eqn2 (Excel)
% Average of Kp1-3 to fit more in the PID characteristic relationship
Ki = max((Ki1 + Ki2 + Ki3) / 3, 0);

% 2.C. Generation of Step Response Graph
num=[6.68*Kd 6.68*Kp 6.68*Ki];
den=[114.89 91.44+(6.68*Kd) 32+(6.68*Kp) (6.68*Ki)];
t=0:0.01:5;
step(num,den,t);
sys = tf(num,den);
s1 = [stepinfo(sys).RiseTime stepinfo(sys).Overshoot];
s2 = [stepinfo(sys).SettlingTime 1/(0.20875*Ki)];
fprintf('s = Rise Time\tPercent Overshoot\tSettling Time\tSteady-State Error')
s = [s1 s2]

% 3. Routh Table
syms EPS
ra=routh(den, EPS)

% 4. Pole Zero Map
%{
sys = tf(num,den);
h = pzplot(sys);
%}

% 5. Transfer Function to State Space Representation
b = num;
a = den;
[A,B,C,D] = tf2ss(b,a)