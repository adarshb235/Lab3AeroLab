clc;
clear;
close all;


alpha = (-10:15) .* (pi/180);

tipCode = 0012;
panelCount = 64;
tipVec = digi(tipCode);
[x_b_tip, y_b_tip] = genNACA(1, tipVec(1), tipVec(2), str2double(strcat(num2str(tipVec(3)), num2str(tipVec(4)))), (panelCount/2) + 2, 0);

rootCode = 2412;
panelCount_2 = 68;
rootVec = digi(rootCode);
[x_b_root, y_b_root] = genNACA(1, rootVec(1), rootVec(2), str2double(strcat(num2str(rootVec(3)), num2str(rootVec(4)))), (panelCount_2/2) + 2, 0);

alpha_2D_deg = -4:2:4; % Degrees for 2D analysis
alpha_2D_rad = alpha_2D_deg .* (pi/180); 
Cl_tip_2D = zeros(1, length(alpha_2D_deg));
Cl_root_2D = zeros(1, length(alpha_2D_deg));

for i = 1:length(alpha_2D_deg)
    Cl_tip_2D(i)  = Vortex_Panel(x_b_tip, y_b_tip, alpha_2D_deg(i));
    Cl_root_2D(i) = Vortex_Panel(x_b_root, y_b_root, alpha_2D_deg(i));
end

p_tip = polyfit(alpha_2D_rad, Cl_tip_2D, 1);
a0_t = p_tip(1);           % Lift slope [per rad]
aero_t = -p_tip(2)/p_tip(1); % Zero-lift AoA [rad] (x-intercept)

p_root = polyfit(alpha_2D_rad, Cl_root_2D, 1);
a0_r = p_root(1);          % Lift slope [per rad]
aero_r = -p_root(2)/p_root(1); % Zero-lift AoA [rad] (x-intercept)

%wing characteristics
b= 36; %span [feet]
c_r = 5 + 4/12; % chord at root [ft]
c_t = 3 + 7/12; % chord at tips [ft]


N = 5; % # of odd terms for circulation

A_mat = zeros(length(alpha), N);
liftingLine = struct();
Cl_ref = struct();

alpha_exp = [-4 -2 0 2 4 6 8] .* pi/180;
p_drag = [0.00232 0.00144 0.00114 0.00144 0.00232 0.00395 0.00612]; % fit data from airfoiltools for NACA 0012


for i = 1:length(alpha)
    alpha_1 = alpha(i);
    geo_r = alpha(i) + (2 * pi / 180); % geometric AoA (geometric twist + alpha) at tips [radians]
    geo_t = alpha(i); % geometric AoA (geometric twist + alpha) at root [radians]
    Cl_ref.tip(i) = Vortex_Panel(x_b_tip, y_b_tip, alpha(i));
    Cl_ref.root(i) = Vortex_Panel(x_b_root, y_b_root, alpha(i));
    [liftingLine.e(i),liftingLine.CL(i),liftingLine.CD_i(i), A_mat(i, :)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r, N);
    
    liftFit = polyfit(alpha_exp, p_drag, 2);
    liftingLine.CD_prof(i) = polyval(liftFit, alpha(i));
    liftingLine.CD_tot(i) = liftingLine.CD_prof(i) + liftingLine.CD_i(i);
end
S = (c_r + c_t)/2 * b;


figure(); % Task 1 Plot
hold on;
grid on;
box on;
plot(alpha * 180/pi, liftingLine.CL, 'LineWidth', 1.5);
xlabel('Angle of Attack (deg)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Angle of Attack');
legend('Cessna 180 (NACA 2412 Root, NACA 0012 Tip)', 'Location', 'best')

figure; % task 2
hold on;
grid on;
box on;
plot(alpha * 180/pi, liftingLine.CD_prof, 'LineWidth', 1.5);
plot(alpha_exp * 180/pi, p_drag, 'LineWidth', 1.2)
xlabel('Angle of Attack (deg)');
ylabel('Drag Coefficient (C_D)');
title('Drag Coefficient vs Angle of Attack');
legend('Calc Profile Drag', 'Experimental Data', 'Location', 'best')

figure; % task 3 plot
hold on;
grid on;
box on;
plot(alpha * 180/pi, liftingLine.CD_tot, 'k', 'LineWidth', 2);
plot(alpha * 180/pi, liftingLine.CD_i, 'b--', 'LineWidth', 1.5);
plot(alpha * 180/pi, liftingLine.CD_prof, 'r-.', 'LineWidth', 1.5);
xlabel('Angle of Attack (deg)');
ylabel('Drag Coefficient');
legend('Total Drag (C_D)', 'Induced Drag (C_{Di})', 'Profile Drag (C_{D0})', 'Location', 'best');
title('Drag Components vs Angle of Attack');

rho = 0.001756; % density [slugs/ft^3] at 10,000' standard day
W = 2500; % aircraft weight [lbs]

V = zeros(1, length(alpha));
T_req = zeros(1, length(alpha));

for i = 1:length(alpha)
    if liftingLine.CL(i) > 0
        V_1 = sqrt((2*W)/(rho * S * liftingLine.CL(i)));
        V(i) = 0.592484 * V_1;
        T_req(i) = 0.5 * rho * V_1^2 * S * liftingLine.CD_tot(i);
    else
        V(i) = NaN;
        T_req(i) = NaN;
    end
end

figure(); % task 4 plot
plot(V, T_req, 'b-', 'LineWidth', 2);
xlabel('Airspeed (knots)');
ylabel('Thrust Required (lbs)');
title('Thrust Required for Steady Level Flight (10,000 ft)');
legend('Cessna 180 (NACA 2412 Root, NACA 0012 Tip)', 'Location', 'best')
grid on;

filenamesVec = ["task1liftpolar", "task2dragvsangle", "task3dragpolar", "task4velvsthrust"];

% Part 1 NACA Code Generator

function [x_b, y_b] = genNACA(c, m, p, t, panel, graph)

m = m/100; 
p = p/10;
t = t/100;

theta = linspace(0, pi, panel);

x = c./2 * (1-cos(theta));
x_c = x/c;

y_t = (t/0.2) * c * ((0.2969 * sqrt(x_c)) - (0.126 * x_c) - (0.3516 * (x_c).^2) + (0.2843 * (x_c).^3) - (0.1036 * (x_c).^4));

y_c = zeros(1, length(x));

if (m == 0 || p == 0)
    dyc_dx = 0;
else
    for i = 1:length(x)
        if (x(i) >= 0) && (x(i) < p*c)
            y_c(i) = m * (x(i) / p^2) * (2*p - x_c(i));
        elseif (x(i) >= p*c) && (x(i) <= c)
            y_c(i) = m * ((c - x(i)) / (1 - p)^2) * (1 + (x_c(i)) - 2*p);
        end
    end
    dyc_dx = gradient(y_c, x);
end


xi = atan(dyc_dx);

x_u = x - y_t .* sin(xi);
x_l = x + y_t .* sin(xi);
y_u = y_c + y_t .* cos(xi);
y_l = y_c - y_t .* cos(xi);

x_b = [fliplr(x_l), x_u(2:end)];
y_b = [fliplr(y_l), y_u(2:end)];

if (graph == 1)
    figure;
    hold on;
    axis equal;
    grid on;
    box on;
    plot(x_b, y_b, 'k');
    plot(x, y_c)
    title("NACA" + m*100 + p*10 + t*100)
    legend('Foil Outline', 'Mean Camber Line', 'Location', 'best');
    filename = append('airfoilNACA', num2str(m*100), num2str(p*10), num2str(t*100)); 
    ax = gca; 
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    print(filename,'-r500','-dpng')
end

end

% Part 2 PLLT Code

function [e,c_L,c_Di,A_matrix] = PLLT(b, a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    theta = (1:N)' * pi / (2*N); % setup and treating data
    c = c_r + (c_t-c_r) .* cos(theta);
    a0 = a0_r + (a0_t - a0_r) .* cos(theta);
    alphaL0 = aero_r + (aero_t - aero_r) .* cos(theta);
    alpha_geo = geo_r + (geo_t - geo_r) .* cos(theta);
    b_i = (alpha_geo - alphaL0);
    Mij = zeros(N);
    for i = 1:N % iterate
        ac_Theta = theta(i);
        mu = 4 * b / (a0(i) * c(i));
        for j = 1:N % iterate
            n = 2 * j - 1;
            Mij(i,j) = (mu + (n / sin(ac_Theta))) * sin(n * ac_Theta);
        end
    end
    A_matrix = Mij \ b_i; % matrix out
    
    S = b*((c_r +c_t)/2); % geometric setup
    AR = b^2/S;
    c_L = A_matrix(1)*pi*AR; %lift coeff

    sum_n_An2 = 0; % iterative solve
    for j=1:N
        n = 2 * j - 1;
        sum_n_An2 = sum_n_An2 + n * (A_matrix(j).^2);
    end

   c_Di = pi * AR * sum_n_An2; % induced drag

   if sum_n_An2 == 0
        e = 0;
    else
        e = A_matrix(1)^2 / sum_n_An2;
    end
    
    
end

% Aux Function To Supplement Use of Part 1 Work

function [vec] = digi(input)
    vec = num2str(input);
    if length(vec) < 4
        % Pad with leading zeros if it's shorter than 4
        vec = pad(vec, 4, 'left', '0'); 
    elseif length(vec) > 4
        % Truncate to 4 characters if it's longer
        vec = vec(1:4);
    end
    vec = vec - '0';
end
