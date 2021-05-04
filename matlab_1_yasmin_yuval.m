%these 3 lines for clearing and closing all running plots
clc;                               
close all;
clear all;
 
%ID of Yuval for the parameters of the project
A=3;
B=1;
C=2;
D=2;
E=3;
F=8;
G=2;
H=0;
I=7;
 
%variables for the project
freq = 10^9;
w = 2*pi*freq;
Z_c = 50;
v_p = 0.8*3*10^8;
L = 0.75;
t_p = L/v_p; 
z = -L/2;
R_g = 20+E*3^(-B);
R_l = 160+F*2^(-D);
gama_l = (R_l-Z_c)/(R_l+Z_c);
gama_g = (R_g-Z_c)/(R_g+Z_c);
t_range = 0:t_p/500:10*t_p;
Vol_val = 0:t_p/500:10*t_p;
Cur_val = 0:t_p/500:10*t_p;
beta = w/v_p;
lambda = 2*pi/beta;
 
v_forward  = 0:t_p/500:10*t_p;
v_backwards = 0:t_p/500:10*t_p;
imp_ratio = Z_c/(Z_c+R_g);
 
%zeroing voltage and current after each iteration
for n=0:t_p/500:10*t_p
    Vol_val = 0;
    Cur_val = 0;
    v_forward  = 0;
    v_backwards = 0;
end

%% section C - 2D plot of voltage on the line
 
z_range= -L:L/500:0;            %res for z axis
t_range_2 = 0:t_p/50:10*t_p;    %res for z axis

[T_vector,Z_vector]= meshgrid(t_range_2,z_range);
v_forward = zeros(length(z_range),length(t_range_2));
v_backwards = zeros(length(z_range),length(t_range_2));
Vol_val_z_t = zeros(length(z_range),length(t_range_2));
 
for n=1:5 
    %plotting forward wave
    t_plus= T_vector - Z_vector/v_p -L/v_p-2*(n-1)*t_p;
    V_g_plus=sin(w.*t_plus).*heaviside(t_plus);
    v_forward=(gama_l.^(n-1)).*(gama_g.^(n-1)).*imp_ratio.*V_g_plus;
    
    %plotting returning wave
    t_minus=T_vector + Z_vector/v_p -L/v_p-2*(n-1)*t_p;
    V_g_minus=sin(w.*t_minus).*heaviside(t_minus);
    v_backwards=(gama_l.^(n)).*(gama_g.^(n-1)).*imp_ratio.*V_g_minus;
    
    %the voltage is superposition of the two waves
    Vol_val_z_t = Vol_val_z_t + v_forward + v_backwards;
end
 
figure(1)
mesh(T_vector,Z_vector,Vol_val_z_t);
colorbar;
grid on;
xlabel('Time [Sec]');
ylabel('Z [m]');
zlabel('Voltage [v]');
title('V(t,z) as $0 < t < 10{\tau_p}$','interpreter','latex');


%% section E

t = 0:t_p/500:10*t_p; % Time vector
V_vector = zeros(size(t)); % Voltage resonace vector
I_vector = zeros(size(t)); % Current resonace vector
Vs = V_vector;
Is = I_vector;

% Calculating the Voltage and Current vector
for i = 0:5
    % maybe change some names
   V_p = C.*(gama_l.^(i)).*(gama_g.^(i)).*sin(w.*(t-(z/v_p)-(2*(i)+1)*t_p)).*heaviside(t-(z/v_p)-(2*(i)+1)*t_p); % for the V+ vector
   V_m = C.*(gama_l.^(i+1)).*(gama_g.^(i)).*sin(w.*(t+(z/v_p)-(2*(i)+1)*t_p)).*heaviside(t+(z/v_p)-(2*(i)+1)*t_p); % for the V- vector
   V_vector = V_vector+V_p+V_m; % Summing the voltage vector
   I_vector = I_vector +V_p/(Z_c)+V_m/(-Z_c);% Summing the Current vector
end

% Calculating the Voltage and Current vector for steady state
% need to make sure plos/ minus
V_ps = 0.509*cos(2.86*pi+w*t-beta*z); % for the V+ vector
V_ms = 0.268*cos(2.86*pi+w*t+beta*z); % for the V- vector
Vs = Vs+V_ps+V_ms; % Summing the voltage vector
Is = Is+V_ps/(Z_c)+V_ms/(-Z_c);% Summing the Current vector

% Plot Voltage at Z=L/2 0<t<10t_p for Excited and Steady state
figure(2)
subplot(2,1,1);
plot(t,V_vector);
ylabel('Voltage [V]');
xlabel('Time [Sec]');
title('Voltage at $Z=\frac{L}{2}$ for $0<t<10{\tau_p}$ - Excited state','interpreter','latex')
subplot(2,1,2);
plot(t,Vs);
ylabel('Voltage [V]');
xlabel('Time [Sec]');
title('Voltage at $Z=\frac{L}{2}$ for $0<t<10{\tau_p}$ - Steady state','interpreter','latex')

% Plot Current at Z=L/2 0<t<10t_p for Excited and Steady state
figure(3)
subplot(2,1,1);
plot(t,I_vector);
ylabel('Current [A]');
xlabel('Time [Sec]');
title('Current at $Z=\frac{L}{2}$ for $0<t<10{\tau_p}$ - Excited state','interpreter','latex')
subplot(2,1,2);
plot(t,Is);
ylabel('Current [A]');
xlabel('Time [Sec]');
title('Current at $Z=\frac{L}{2}$ for $0<t<10{\tau_p}$ - Steady state','interpreter','latex')

%% Section F

T = 1/freq;
T_vec = [0, 0.1*T, 0.3*T, 0.5*T, 0.6*T, 0.8*T];
Z_axis = -L:L/500:0;

figure(4)
for i=1:6
    V_ps = 0.509*cos(2.86*pi+w*T_vec(i)-beta*Z_axis); % for the V+ vector
    V_ms = 0.268*cos(2.86*pi+w*T_vec(i)+beta*Z_axis); % for the V- vector
    Vs = V_ps+V_ms; % Summing the voltage vector
    plot(Z_axis,Vs);
    hold on
end
ylabel('Voltage [V] ');
xlabel('Z [m]');
title('V(z) for different times');
legend('t=0T','t=0.1T','t=0.3T','t=0.5T','t=0.6T','t=0.8T')
hold off;

%% bonus

t_points = 7.5*30; 
z_points = 1000;
z = linspace(-L, 0, z_points);
t = linspace(0, 10*t_p, t_points)';
W = [];
for i = z
    [v_right,v_left] = wave_init(t,i,v_p, L, R_g, w, Z_c, gama_l, gama_g);
    W = [W , v_right + v_left];
end

figure(5);
for i = [1:1:t_points]
    plot(z,W(i,:));
    title("V(t) on the line at t="+t(i)+" sec");
    xlim([-0.75,0]);
    xlabel("z[m]");
    ylim([-1.5,1.5]);
    ylabel("voltage[V]");
    pause(1/5);
end


function [v_right,v_left] = wave_init(t, z, v_p, L, R_g, w, Z_c, gama_l, gama_g)

    step = @(x) x>=0; 
    %global vp; global Rg; global Rl; global l; global omega; global Zc; 
    t_p = L/v_p;
    %gama_g = (Rg-Zc)/(Rg+Zc);
    %gama_l = (Rl-Zc)/(Rl+Zc);
    v_right = zeros(size(t,1),size(t,2));
    v_left = zeros(size(t,1),size(t,2));
    amp = Z_c/(Z_c+R_g); %a factor for reducing the wave amplitude
    
    % add first 2*5 waves on the line (fsource wave an 9 echos)
    for n = [1:1:5]
        arg = t - z/v_p - t_p*(2*n-1);
        v_right = v_right + amp*step(arg).*sin(w*arg);
        amp = amp * gama_l; %add a bounce effect on load side
        
        arg = t + z/v_p - t_p*(2*n-1);
        v_left = v_left + amp*step(arg).*sin(w*arg);
        amp = amp * gama_g; %add a bounce effect on source side
    end
end
